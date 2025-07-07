
def extracting_value_from_string(string, val_name):
    string = string.replace(val_name, '').strip()
    def is_part_of_number(charachter):
        if charachter.isdigit() or charachter == '.':
            return True
        else:
            return False
    return float(''.join(filter(is_part_of_number, string)))



def extracting_values_from_txt_file(file_path, val_names):
    if isinstance(val_names, str):
        val_names = [val_names]
    file = open(file_path, "r")
    data = {}

    for line in file:
        for val_name in val_names:
            if line.startswith(val_name):
                if val_name in data:
                    raise AssertionError(val_name + ' is listed more than once in file.')
                data[val_name] = extracting_value_from_string(line, val_name)
    
    return data
