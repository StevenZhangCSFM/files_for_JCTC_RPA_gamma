import re
from collections import defaultdict

# Function to parse the log file and sum the times for each category
def parse_and_sum_times(file_path):
    # Dictionary to store the summed times for each part
    time_sums = defaultdict(float)
    
    # Regular expression pattern to match the different parts and their times
    pattern = re.compile(r"sumOperatorMulTime (\d+\.\d+) ms, sumYT_mul_YTime (\d+\.\d+) ms, sum_orthonormalTime (\d+\.\d+) ms, sumYT_operator_YTime (\d+\.\d+) ms, sumTransfer_HpTime (\d+\.\d+) ms, sumSpectral_quadratureTime (\d+\.\d+) ms")
    
    with open(file_path, 'r') as file:
        for line in file:
            match = pattern.search(line)
            if match:
                time_sums["sumOperatorMulTime"] += float(match.group(1))
                time_sums["sumYT_mul_YTime"] += float(match.group(2))
                time_sums["sum_orthonormalTime"] += float(match.group(3))
                time_sums["sumYT_operator_YTime"] += float(match.group(4))
                time_sums["sumTransfer_HpTime"] += float(match.group(5))
                time_sums["sumSpectral_quadratureTime"] += float(match.group(6))
    
    return time_sums

# Function to display the total times
def display_total_times(time_sums):
    print(f"In total, sumOperatorMulTime {time_sums['sumOperatorMulTime']:.1f} ms, "
          f"sumYT_mul_YTime {time_sums['sumYT_mul_YTime']:.1f} ms, "
          f"sum_orthonormalTime {time_sums['sum_orthonormalTime']:.1f} ms, "
          f"sumYT_operator_YTime {time_sums['sumYT_operator_YTime']:.1f} ms, "
          f"sumTransfer_HpTime {time_sums['sumTransfer_HpTime']:.1f} ms, "
          f"sumSpectral_quadratureTime {time_sums['sumSpectral_quadratureTime']:.1f} ms")

# Main function to run the analysis
def main():
    file_path = "out.log"  # Replace with the path to your input file
    time_sums = parse_and_sum_times(file_path)
    display_total_times(time_sums)

if __name__ == "__main__":
    main()