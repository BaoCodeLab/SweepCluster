from library import parameters

if __name__ == "__main__":
    args=parameters.parse_arguments()
    args.func(args)
