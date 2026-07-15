from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib.pyplot as plt
from obspy import read as obspy_read


def load_geometry_file(path: Path) -> np.ndarray:
	positions: list[float] = []
	with open(path, encoding="utf-8") as fh:
		for line in fh:
			s = line.strip()
			if not s or s.lower().startswith("trace"):
				continue
			parts = s.replace(",", ".").split()
			if len(parts) < 2:
				continue
			try:
				positions.append(float(parts[1]))
			except ValueError:
				continue
	if not positions:
		raise ValueError(f"No receiver positions read from geometry file: {path}")
	return np.asarray(positions, dtype=float)


def ormsby_response(freqs: np.ndarray, f1: float, f2: float,
					f3: float, f4: float) -> np.ndarray:
	h = np.zeros_like(freqs, dtype=float)
	m_lo = (freqs > f1) & (freqs <= f2)
	m_pb = (freqs > f2) & (freqs <= f3)
	m_hi = (freqs > f3) & (freqs < f4)
	h[m_lo] = (freqs[m_lo] - f1) / max(f2 - f1, 1e-12)
	h[m_pb] = 1.0
	h[m_hi] = (f4 - freqs[m_hi]) / max(f4 - f3, 1e-12)
	return h


def ormsby(trace: np.ndarray, dt_s: float, f1: float, f2: float,
		   f3: float, f4: float, pad_pct: float = 0.25) -> np.ndarray:
	n = len(trace)
	n_pad = int(n * (1.0 + pad_pct))
	n_fft = 1 << (max(n_pad, 2) - 1).bit_length()
	spec = np.fft.rfft(trace.astype(np.float64), n=n_fft)
	freqs = np.fft.rfftfreq(n_fft, d=dt_s)
	h = ormsby_response(freqs, f1, f2, f3, f4)
	out = np.fft.irfft(spec * h, n=n_fft)
	return out[:n].astype(np.float32)


def apply_ormsby_all(data: np.ndarray, dt_s: float,
					 f1: float, f2: float, f3: float, f4: float,
					 pad_pct: float = 0.25) -> np.ndarray:
	out = np.empty_like(data)
	for i in range(data.shape[0]):
		out[i] = ormsby(data[i], dt_s, f1, f2, f3, f4, pad_pct=pad_pct)
	return out


def read_seg2(seg2_path: Path) -> dict[str, Any]:
	st = obspy_read(str(seg2_path), format="SEG2")
	n_traces = len(st)
	n_samp = max(tr.stats.npts for tr in st)
	dt_s = float(st[0].stats.delta)

	data = np.zeros((n_traces, n_samp), dtype=np.float32)
	hdrs: list[dict[str, Any]] = []

	recv_locs: list[float | None] = []
	for i, tr in enumerate(st):
		npts = tr.stats.npts
		data[i, :npts] = tr.data.astype(np.float32)
		seg2_hdr = dict(getattr(tr.stats, "seg2", {}))
		hdrs.append(seg2_hdr)

		loc = str(seg2_hdr.get("RECEIVER_LOCATION", "")).strip().replace(",", ".")
		recv_locs.append(float(loc) if loc else None)

	hdr0 = hdrs[0] if hdrs else {}
	delay_s = float(str(hdr0.get("DELAY", "0")).replace(",", "."))
	delay_ms = delay_s * 1000.0

	shot_pos = None
	src = str(hdr0.get("SOURCE_LOCATION", "")).strip().replace(",", ".")
	if src:
		try:
			shot_pos = float(src.split()[0])
		except ValueError:
			shot_pos = None

	ffid = None
	ssn = str(hdr0.get("SHOT_SEQUENCE_NUMBER", "")).strip()
	if ssn.isdigit():
		ffid = int(ssn)

	recv_locs_arr: np.ndarray | None = None
	if all(v is not None for v in recv_locs):
		recv_locs_arr = np.asarray(recv_locs, dtype=float)

	time_ms = delay_ms + np.arange(n_samp) * dt_s * 1000.0

	return {
		"data": data,
		"dt_s": dt_s,
		"n_traces": n_traces,
		"n_samp": n_samp,
		"delay_ms": delay_ms,
		"shot_pos_m": shot_pos,
		"ffid": ffid,
		"recv_locs_m": recv_locs_arr,
		"time_ms": time_ms,
		"trace_headers": hdrs,
		"global_header": hdr0,
	}


def print_header_summary(seg2_path: Path, info: dict[str, Any],
						 max_trace_headers: int = 5) -> None:
	print("\n=== SEG2 SUMMARY ===")
	print(f"File             : {seg2_path}")
	print(f"FFID             : {info['ffid']}")
	print(f"Traces           : {info['n_traces']}")
	print(f"Samples/trace    : {info['n_samp']}")
	print(f"dt               : {info['dt_s'] * 1000.0:.6f} ms")
	print(f"DELAY            : {info['delay_ms']:.3f} ms")
	print(f"SOURCE_LOCATION  : {info['shot_pos_m']}")

	gh = info["global_header"]
	keys = [
		"ACQUISITION_DATE", "ACQUISITION_TIME", "CLIENT", "COMPANY",
		"LINE_ID", "SHOT_SEQUENCE_NUMBER", "DELAY", "SOURCE_LOCATION",
		"SAMPLE_INTERVAL", "STACK", "INSTRUMENT",
	]
	print("\n--- Selected global header fields ---")
	for k in keys:
		if k in gh:
			print(f"{k:20s}: {gh.get(k)}")

	print("\n--- First trace header snippets ---")
	n_show = min(max_trace_headers, len(info["trace_headers"]))
	for i in range(n_show):
		h = info["trace_headers"][i]
		print(
			f"Trace {i+1:02d}: "
			f"CHAN={h.get('CHANNEL_NUMBER', 'n/a')}  "
			f"RECV_LOC={h.get('RECEIVER_LOCATION', 'n/a')}  "
			f"TRACE_ID={h.get('TRACE_ID', 'n/a')}"
		)


def save_data_exports(out_dir: Path, base_name: str,
					  info: dict[str, Any], data_filt: np.ndarray,
					  offsets_m: np.ndarray) -> None:
	out_dir.mkdir(parents=True, exist_ok=True)
	time_ms = info["time_ms"]
	data = info["data"]

	txt_raw = out_dir / f"{base_name}_raw.txt"
	txt_filt = out_dir / f"{base_name}_filt.txt"
	dat_raw = out_dir / f"{base_name}_raw.dat"
	dat_filt = out_dir / f"{base_name}_filt.dat"
	meta_json = out_dir / f"{base_name}_meta.json"

	hdr_cols = ["time_ms"] + [f"tr{i+1}" for i in range(data.shape[0])]
	arr_raw = np.column_stack([time_ms, data.T])
	arr_filt = np.column_stack([time_ms, data_filt.T])

	np.savetxt(txt_raw, arr_raw, fmt="%.6f", header=" ".join(hdr_cols), comments="")
	np.savetxt(txt_filt, arr_filt, fmt="%.6f", header=" ".join(hdr_cols), comments="")

	data.astype(np.float32).tofile(dat_raw)
	data_filt.astype(np.float32).tofile(dat_filt)

	meta = {
		"n_traces": int(info["n_traces"]),
		"n_samp": int(info["n_samp"]),
		"dt_s": float(info["dt_s"]),
		"delay_ms": float(info["delay_ms"]),
		"shot_pos_m": info["shot_pos_m"],
		"offsets_m": offsets_m.tolist(),
		"dat_raw": dat_raw.name,
		"dat_filt": dat_filt.name,
		"txt_raw": txt_raw.name,
		"txt_filt": txt_filt.name,
	}
	with open(meta_json, "w", encoding="utf-8") as fh:
		json.dump(meta, fh, indent=2)

	print("\n=== Saved data ===")
	print(f"TXT raw          : {txt_raw}")
	print(f"TXT filtered     : {txt_filt}")
	print(f"DAT raw (f32)    : {dat_raw}")
	print(f"DAT filtered(f32): {dat_filt}")
	print(f"Meta JSON        : {meta_json}")


def _norm_for_wiggles(data: np.ndarray) -> np.ndarray:
	stds = np.asarray([data[i].astype(float).std() for i in range(data.shape[0])], dtype=float)
	valid = stds[stds > 1e-20]
	med = float(np.median(valid)) if len(valid) else 1.0
	norms = np.clip(stds, med * 0.3, med * 3.0)
	norms = np.where(norms > 1e-20, norms, med)
	return norms


def apply_gain(data: np.ndarray, dt_s: float, mode: str,
			   window_ms: float, stat: str) -> np.ndarray:
	mode_l = mode.lower()
	if mode_l == "none":
		return data.copy()

	out = data.astype(np.float32).copy()

	if mode_l == "norm":
		for i in range(out.shape[0]):
			mx = float(np.max(np.abs(out[i])))
			if mx > 1e-20:
				out[i] /= mx
		return out

	if mode_l == "agc":
		win = max(3, int((window_ms / 1000.0) / dt_s))
		if win % 2 == 0:
			win += 1
		ker = np.ones(win, dtype=np.float64) / float(win)

		for i in range(out.shape[0]):
			tr = out[i].astype(np.float64)
			if stat.lower() == "mean":
				env = np.convolve(np.abs(tr), ker, mode="same")
			else:
				env = np.sqrt(np.convolve(tr * tr, ker, mode="same"))
			env = np.where(env > 1e-12, env, 1e-12)
			out[i] = (tr / env).astype(np.float32)
		return out

	raise ValueError(f"Unknown gain mode: {mode}")


def estimate_three_shots(recv_x: np.ndarray) -> dict[str, float]:
	dx_start = float(recv_x[1] - recv_x[0]) if len(recv_x) > 1 else 1.0
	dx_end = float(recv_x[-1] - recv_x[-2]) if len(recv_x) > 1 else 1.0
	mid_lo = float(recv_x[len(recv_x) // 2 - 1])
	mid_hi = float(recv_x[len(recv_x) // 2])
	return {
		"first": float(recv_x[0] - dx_start),
		"middle": float((mid_lo + mid_hi) / 2.0),
		"last": float(recv_x[-1] + dx_end),
	}


def resolve_shot_position(shot_pos_header: float | None,
				  recv_x: np.ndarray,
				  file_idx_1based: int,
				  total_files: int,
				  mode: str) -> tuple[float, str]:
	mode_l = mode.lower()
	if mode_l == "header" and shot_pos_header is not None:
		return float(shot_pos_header), "header"

	est = estimate_three_shots(recv_x)
	if total_files <= 1:
		return est["middle"], "estimated-middle"

	if file_idx_1based == 1:
		return est["first"], "estimated-first"
	if file_idx_1based == total_files:
		return est["last"], "estimated-last"
	if file_idx_1based == (total_files // 2 + (1 if total_files % 2 else 0)):
		return est["middle"], "estimated-middle"

	if shot_pos_header is not None:
		return float(shot_pos_header), "header-fallback"
	return est["middle"], "estimated-middle-fallback"


def plot_against_offset(out_dir: Path, base_name: str,
						info: dict[str, Any], data_filt: np.ndarray,
						offsets_m: np.ndarray,
						dt_s: float,
						gain_mode: str,
						agc_window_ms: float,
						agc_stat: str,
						show_plots: bool) -> None:
	out_dir.mkdir(parents=True, exist_ok=True)

	data_raw = info["data"]
	time_ms = info["time_ms"]
	time_no_delay_ms = time_ms - float(info["delay_ms"])

	raw_gain = apply_gain(data_raw, dt_s, gain_mode, agc_window_ms, agc_stat)
	fil_gain = apply_gain(data_filt, dt_s, gain_mode, agc_window_ms, agc_stat)

	x = offsets_m
	x_margin = (float(np.ptp(x)) / max(len(x) - 1, 1)) * 0.6 if len(x) > 1 else 1.0

	def _plot_wiggles(t_axis: np.ndarray, suffix: str, axis_label: str) -> Path:
		fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
		for ax, data, title in [
			(axes[0], raw_gain, "Raw"),
			(axes[1], fil_gain, "Filtered (Ormsby)"),
		]:
			norms = _norm_for_wiggles(data)
			dx = float(np.ptp(x)) / max(len(x) - 1, 1) if len(x) > 1 else 1.0
			scale = norms * 2.0
			for i, off in enumerate(x):
				trn = np.clip(data[i].astype(float) / scale[i], -1.0, 1.0)
				ax.plot(off + trn * dx * 0.85, t_axis, color="black", lw=0.5, alpha=0.85)
			ax.axhline(0.0, color="green", ls=":", lw=0.8, alpha=0.8)
			ax.set_title(f"{title} | gain={gain_mode}")
			ax.set_xlabel("Offset from shot (m)")
			ax.set_xlim(x.min() - x_margin, x.max() + x_margin)
			ax.set_ylim(float(np.max(t_axis)), float(np.min(t_axis)))
			ax.grid(True, lw=0.3, alpha=0.3)

		axes[0].set_ylabel(axis_label)
		wiggle_png = out_dir / f"{base_name}_offset_wiggles_{suffix}.png"
		fig.suptitle(f"{base_name} -- offsets/wiggles ({suffix})")
		fig.savefig(wiggle_png, dpi=180)
		if show_plots:
			fig.show()
		else:
			plt.close(fig)
		return wiggle_png

	def _plot_image(t_axis: np.ndarray, suffix: str, axis_label: str) -> Path:
		fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
		vmax = np.percentile(np.abs(fil_gain), 98)
		img = ax.imshow(
			fil_gain.T,
			aspect="auto",
			cmap="seismic",
			vmin=-vmax,
			vmax=vmax,
			extent=[x.min(), x.max(), t_axis[-1], t_axis[0]],
			interpolation="nearest",
		)
		ax.axhline(0.0, color="k", ls=":", lw=0.8)
		ax.set_xlabel("Offset from shot (m)")
		ax.set_ylabel(axis_label)
		ax.set_title(f"Filtered gather (image view) | gain={gain_mode} | {suffix}")
		fig.colorbar(img, ax=ax, label="Amplitude")
		img_png = out_dir / f"{base_name}_offset_image_{suffix}.png"
		fig.savefig(img_png, dpi=180)
		if show_plots:
			fig.show()
		else:
			plt.close(fig)
		return img_png

	wiggle_with = _plot_wiggles(time_ms, "with_delay", "Time (ms, with delay)")
	image_with = _plot_image(time_ms, "with_delay", "Time (ms, with delay)")
	wiggle_no = _plot_wiggles(time_no_delay_ms, "no_delay", "Time (ms, no delay)")
	image_no = _plot_image(time_no_delay_ms, "no_delay", "Time (ms, no delay)")

	if show_plots:
		plt.show()

	print("\n=== Saved plots ===")
	print(f"Wiggles PNG      : {wiggle_with}")
	print(f"Image PNG        : {image_with}")
	print(f"Wiggles PNG      : {wiggle_no}")
	print(f"Image PNG        : {image_no}")


def _ffid_from_name(path: Path) -> int:
	digits = "".join(ch for ch in path.stem if ch.isdigit())
	return int(digits) if digits else 0


def _collect_seg2_files(input_path: Path) -> list[Path]:
	if input_path.is_file():
		return [input_path]
	if not input_path.is_dir():
		return []
	candidates = list(input_path.glob("*.seg2")) + list(input_path.glob("*.SEG2"))
	unique = {str(p.resolve()).lower(): p for p in candidates}
	return sorted(unique.values(), key=_ffid_from_name)


def _default_out_dir(input_path: Path, out_name: str | None) -> Path:
	# Place outputs in outer lvl/output as requested.
	# Example input folder .../lvl/data/150 -> .../lvl/output/150_
	if input_path.is_dir():
		base = input_path.name
		lvl_root = input_path.parent.parent
	elif input_path.is_file():
		base = input_path.parent.name
		lvl_root = input_path.parent.parent
	else:
		base = "readseg"
		lvl_root = Path.cwd()
	folder_name = out_name if out_name else f"{base}_"
	return lvl_root / "output" / folder_name


def process_one_file(seg2_path: Path, out_dir: Path,
			 file_idx_1based: int, total_files: int,
			 args: argparse.Namespace) -> None:
	base_name = seg2_path.stem

	info = read_seg2(seg2_path)
	print_header_summary(seg2_path, info, max_trace_headers=args.max_trace_headers)

	recv_x: np.ndarray
	if args.geometry_file:
		geom = load_geometry_file(Path(args.geometry_file))
		recv_x = geom[:info["n_traces"]]
		print(f"\nUsing geometry file receiver x: {args.geometry_file}")
	elif info["recv_locs_m"] is not None:
		recv_x = info["recv_locs_m"]
		print("\nUsing SEG2 RECEIVER_LOCATION for receiver x.")
	else:
		recv_x = np.arange(info["n_traces"], dtype=float)
		print("\nNo receiver x in SEG2; using trace index as x.")

	shot_pos, shot_pos_src = resolve_shot_position(
		info["shot_pos_m"], recv_x, file_idx_1based, total_files, args.shot_pos_mode
	)
	print(f"\nShot position source: {shot_pos_src}")

	data_filt = apply_ormsby_all(
		info["data"],
		info["dt_s"],
		args.f1,
		args.f2,
		args.f3,
		args.f4,
		pad_pct=args.pad_pct,
	)

	offsets_m = recv_x - float(shot_pos)

	print("\n=== Offset summary ===")
	print(f"Shot position (m): {shot_pos:.3f}")
	print(f"Receiver x min/max: {recv_x.min():.3f} / {recv_x.max():.3f}")
	print(f"Offset min/max (m): {offsets_m.min():.3f} / {offsets_m.max():.3f}")

	save_data_exports(out_dir, base_name, info, data_filt, offsets_m)
	plot_against_offset(
		out_dir,
		base_name,
		info,
		data_filt,
		offsets_m,
		info["dt_s"],
		args.gain,
		args.agc_window_ms,
		args.agc_stat,
		args.show_plots,
	)


def main() -> None:
	parser = argparse.ArgumentParser(
		description=(
			"Read SEG2 file(s), print key headers, export raw/filtered traces "
			"as txt/dat, and plot data against offset."
		)
	)
	parser.add_argument(
		"input_path",
		help=(
			"Path to one .seg2 file/folder OR profile name under data/ "
			"(e.g. 120)"
		),
	)
	parser.add_argument(
		"geometry",
		nargs="?",
		default=None,
		help="Optional geometry shorthand for profile mode: 100 or 200",
	)
	parser.add_argument("--out-dir", default=None,
						help="Output folder (default: outer lvl/output/<folder>_)")
	parser.add_argument("--out-name", default=None,
						help="Output subfolder name under outer output/ (e.g. 150_)")
	parser.add_argument("--geometry-file", default=None,
						help="Optional geometry file (e.g., geometry100.txt) to override receiver x")
	parser.add_argument("--f1", type=float, default=2.0)
	parser.add_argument("--f2", type=float, default=4.0)
	parser.add_argument("--f3", type=float, default=140.0)
	parser.add_argument("--f4", type=float, default=180.0)
	parser.add_argument("--pad-pct", type=float, default=0.25)
	parser.add_argument("--shot-pos-mode", choices=["estimated", "header"],
						default="estimated",
						help="Shot position source: estimated first/mid/last or SEG2 header")
	parser.add_argument("--gain", choices=["none", "norm", "agc"], default="none",
						help="Gain mode for plotting")
	parser.add_argument("--agc-window-ms", type=float, default=20.0,
						help="AGC window length in ms")
	parser.add_argument("--agc-stat", choices=["mean", "rms"], default="rms",
						help="AGC statistic")
	parser.add_argument("--show-plots", action="store_true", default=True,
						help="Pop up saved figures (with and without delay) [default: on]")
	parser.add_argument("--no-show-plots", action="store_false", dest="show_plots",
						help="Disable popup windows (useful for long batch runs)")
	parser.add_argument("--max-trace-headers", type=int, default=5)
	args = parser.parse_args()

	lvl_root = Path(__file__).resolve().parents[1]
	data_root = lvl_root / "data"

	input_path = Path(args.input_path)
	if not input_path.exists():
		# Compatibility mode: allow `python readseg.py 120 100`
		# where 120 is interpreted as data/<profile>.
		profile_dir = data_root / str(args.input_path)
		if profile_dir.exists():
			input_path = profile_dir
			if args.geometry and not args.geometry_file:
				g = str(args.geometry).strip().lower().replace("geometry", "")
				if g in {"100", "200"}:
					args.geometry_file = str(data_root / f"geometry{g}.txt")
				else:
					raise ValueError(
						f"Invalid geometry shorthand '{args.geometry}'. Use 100 or 200."
					)
		else:
			raise FileNotFoundError(f"Input path/profile not found: {args.input_path}")

	seg2_files = _collect_seg2_files(input_path)
	if not seg2_files:
		raise FileNotFoundError(f"No SEG2 files found at: {input_path}")

	out_dir = Path(args.out_dir) if args.out_dir else _default_out_dir(input_path, args.out_name)
	out_dir.mkdir(parents=True, exist_ok=True)

	print("\n=== Batch run ===")
	print(f"Input            : {input_path}")
	print(f"SEG2 files found : {len(seg2_files)}")
	print(f"Output folder    : {out_dir}")

	for idx, seg2_path in enumerate(seg2_files, start=1):
		print(f"\n{'='*18} FILE {idx}/{len(seg2_files)} {'='*18}")
		process_one_file(seg2_path, out_dir, idx, len(seg2_files), args)

	print("\nDone.")


if __name__ == "__main__":
	main()
