file_name = 'scipy/integrate/_quadpack_py.py'

with open(file_name, 'r') as file:
    for line_number, line in enumerate(file):
        for character_number, character in enumerate(line):
            if character == "\x07":
                print(line_number, character_number)
