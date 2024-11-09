"""Please refer to section 1 of the report and make sure all the libraries are imported"""

# class form beautifulshop library to parse through html files 
from bs4 import BeautifulSoup
# to interact with the operating system.
import os
# finds files and directories by matching patterns 
import glob


# path to the directory with the stored fastqc files that we generated with run_fastqc.sh
path_to_fastqc_files = "./fastqc_results"

path_to_results = "./"

# Output file for results
output_file = os.path.join(path_to_results, "fastqc_ResultCollection.txt")

def parse_fastqc_report(filepath):
    """function to rcollect the metric that are demanded 
    input : path to file (string) """
    with open(filepath, "r") as file:
        # parse the file as an HTML document
        soup = BeautifulSoup(file, "html.parser")

    # initialise directory tot store the collected data points 
    data = {
        # total sequences that are analysrd in the report 
        "total_sequences": None,
        # keep track of quality scores 
        "per_base_quality": [],
        # to check if adapter sequences are present
        "adapter_content": False,
        # check for - and store potential warnings in the file 
        "warnings": []
    }

    # searches in the HTML for an <a> tag that has an href '#M1'
    per_base_quality_section = soup.find("a", href="#M1")
    if per_base_quality_section:
        quality_scores = []
        # located the table we are looking for 
        for row in per_base_quality_section.find_next("table").find_all("tr")[1:]:
            # store the value on the first cell as the position (saved as text)
            position = row.find_all("td")[0].text
            # stores value of the second cell as the mean quality score (saved as text without space)
            mean_quality_text = row.find_all("td")[1].text.strip() 
            
            # conversion into numeric value 
            if mean_quality_text.replace('.', '', 1).isdigit():
                mean_quality = float(mean_quality_text)
                # append position and mean qialuty score to the quality scores list
                quality_scores.append((position, mean_quality))
        # after iteration throught all the rows, the values are added to the data directory under 
        # the 'per_base_quality'-key        
        data["per_base_quality"] = quality_scores

        # Same spiel again this time with #M6 to collect 
        adapter_section = soup.find("a", href="#M6")
        if adapter_section and "FAIL" in adapter_section.find_next("td").text:
            data["adapter_content"] = True

        # save potential warnings under warnings 
        for warning in soup.find_all("li", {"class": "warn"}):
            data["warnings"].append(warning.text)

        return data


# to get all *_fastqc.html files and group them by sample (each two files)
html_files = glob.glob(os.path.join(path_to_fastqc_files, "*_fastqc.html"))
file_pairs = {}

for file_path in html_files:
    # get the name of the file along (both)
    filename = os.path.basename(file_path)
    sample_name, mate = filename.split("_fastqc")[0].rsplit("_", 1)
    
    if sample_name not in file_pairs:
        file_pairs[sample_name] = {}
    
    # assign numbers 1 and 2 to the files respectively 
    file_pairs[sample_name][mate] = file_path

# open the file to store the output in 
with open(output_file, "w") as f:
    # Analyze each sample's mate pairs
    for sample, mates in file_pairs.items():
        if "1" in mates and "2" in mates:  
            data_1 = parse_fastqc_report(mates["1"])
            data_2 = parse_fastqc_report(mates["2"])
            
            # add the total seqquneces to the resultsfile 
            f.write(f"Results for {sample}:\n")
            f.write(f"  Total Sequences - Mate 1: {data_1['total_sequences']}, Mate 2: {data_2['total_sequences']}\n")
            
            # add the adapter content reuslts to the resultsfile 
            adapter_content = data_1["adapter_content"] or data_2["adapter_content"]
            f.write(f"  Adapter Content Detected: {'Yes' if adapter_content else 'No'}\n")
            
            # add thpotential warnings to the resultsfile 
            warnings = data_1["warnings"] + data_2["warnings"]
            f.write(f"  Issues/Warnings: {', '.join(warnings) if warnings else 'None'}\n")
            
            # finally, calculate the average base quality and add them to the results file as well 
            if data_1['per_base_quality'] and data_2['per_base_quality']:
                avg_quality_1 = sum(q[1] for q in data_1["per_base_quality"]) / len(data_1["per_base_quality"])
                avg_quality_2 = sum(q[1] for q in data_2["per_base_quality"]) / len(data_2["per_base_quality"])
                f.write(f"  Average Base Quality - Mate 1: {avg_quality_1:.2f}, Mate 2: {avg_quality_2:.2f}\n")
            
            f.write("\n" + "-"*40 + "\n\n")
        else:
            f.write(f"Warning: Missing mate files for sample {sample}\n\n")

