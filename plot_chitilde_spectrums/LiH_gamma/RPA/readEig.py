import re

def parse_data_file(input_file, output_file):
    # Regular expressions to identify section headers and data rows
    section_header_pattern = re.compile(r"qptIndex (\d+), omegaIndex (\d+)")
    data_row_pattern = re.compile(r"[+-]?\d+\.\d+E[+-]?\d+")
    
    # Dictionary to store data for each section
    sections = {}
    current_section = None
    
    with open(input_file, 'r') as infile:
        for line in infile:
            header_match = section_header_pattern.match(line)
            if header_match:
                # Start a new section
                qpt_index = int(header_match.group(1))
                omega_index = int(header_match.group(2))
                current_section = (qpt_index, omega_index)
                sections[current_section] = []
            elif current_section:
                # Collect data for the current section
                data = data_row_pattern.findall(line)
                if data:
                    sections[current_section].extend(map(float, data))
    
    # Write the output to the file in the desired format
    with open(output_file, 'w') as outfile:
        # Write the header
        header = " ".join([f"({q}, {o})" for q, o in sections.keys()])
        outfile.write(header + "\n")
        
        # Find the maximum number of rows needed
        max_rows = max(len(values) for values in sections.values())
        
        # Write the data row by row
        for i in range(max_rows):
            row = []
            for values in sections.values():
                # Write the value if it exists, otherwise write an empty space
                row.append(f"{values[i]:.6g}" if i < len(values) else "")
            outfile.write(" ".join(row) + "\n")

# File paths
input_file = "LiH.outEig"
output_file = "LiH.txt"

# Parse the data file and create the output
parse_data_file(input_file, output_file)