import argparse


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "mode",
        choices=[
            "desktop",
            "web",
            "refraction"
        ]
    )

    args = parser.parse_args()


    if args.mode == "desktop":
        from apps.desktop.lvl_studio import main
        main()

    elif args.mode == "web":
        import uvicorn
        uvicorn.run(
            "apps.web.server:app",
            host="127.0.0.1",
            port=8000,
            reload=True
        )

    elif args.mode == "refraction":
        from apps.refraction import main
        main()


if __name__ == "__main__":
    main()