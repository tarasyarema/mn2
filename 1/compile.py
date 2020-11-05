def compile():
    file_name = "./YaremaTaras.c"

    input = open(file_name, "r").readlines()
    done = False

    for (i, line) in enumerate(input):
        if done:
            break

        if "#define SUBMIT" in line:
            input[i] = "#define SUBMIT 0"
            done = True

    open("./test.c", "w").writelines(input)


if __name__ == "__main__":
    compile()
