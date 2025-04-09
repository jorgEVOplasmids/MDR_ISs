#!/bin/bash

input_file="$1"
output_file="$2"

temp_file=$(mktemp)

#An awk search, printing first the head of our dataframe. Skip the header for the document and group by the Genome_ID (second column)
awk -F',' '
BEGIN {
    OFS=",";
    print "Taxon ID", "Genome ID", "Genome Name", "Antibiotic", "Resistant Phenotype", "Measurement", 
    "Measurement Sign", "Measurement Value", "Measurement Unit", "Laboratory Typing Method", 
    "Laboratoy Typing Method Version", "Laboratory Typing Platform", "Vendor", "Testing Standard", 
    "Testing Standart Year", "Computational Method", "Computational Method Version", 
    "Computational Method Performance", "Evidence", "Source", "Pubmed" > "'"$temp_file"'"
}
{
    if (NR == 1) next;  
    key = $2;  
    data[key] = data[key] ? data[key] ORS $0 : $0;  
}

END {
    for (key in data) {  
        split(data[key], rows, ORS);  
        
        # Keep the row if there is only one coincidence for Genome_ID and there is information for 
        #the column Antibiotic and Resistance.Phenotype 
        
        if (length(rows) == 1) { 
            split(rows[1], cols, FS);
            if (cols[4] != "" && cols[5] != "" && (cols[5] == "Resistant" || cols[5] == "Susceptible")) {
                print rows[1] >> "'"$temp_file"'";
            }
            
        } else {
            delete Antibiotics;
            delete ValidRows;
            delete Resistant;
            delete Susceptible;
            delete Measurement;
            
            
            min_resistant = 1e10;
            max_susceptible = -1e10;
            min_resistant_row = "";
            min_resistant_row_2 = "";
            max_susceptible_row = "";
            max_susceptible_row_2 = "";
            has_resistant = 0;
            has_susceptible = 0;
            
            # If there are several coincidences for Genome_ID, we extract the minimum value for Resistant
            # M.I.C and the maximum value for Susceptible M.I.C.
            
            for (i in rows) {
                split(rows[i], cols, FS);
            
                if (cols[5] == "Resistant") {
                    has_resistant = 1;
                    valid_resistant = cols[8];
                
                    if (valid_resistant != "") {  
                        if (valid_resistant < min_resistant) {
                            min_resistant = valid_resistant;  
                            min_resistant_row = rows[i];
                        }
                        
                    } else {
                        #min_resistant_row = ""
                        min_resistant_row_2 = rows[i];  
                    }
                    
                } else if (cols[5] == "Susceptible") {
                    has_susceptible = 1;
                    valid_susceptible = cols[8];
                
                    if (valid_susceptible != "") {  
                        if (valid_susceptible > max_susceptible) {
                            max_susceptible = valid_susceptible;  
                            max_susceptible_row = rows[i];
                        }
                    } else {
                        #max_susceptible_row = ""
                        max_susceptible_row_2 = rows[i]; 
                    }
                }
            }
            

            # If there are values both for Resistant and Susceptible data, we keep both if 
            #Resistant > Susceptible data.
            
            if (has_resistant && has_susceptible) {
                if (min_resistant_row != "" && max_susceptible_row != "") {
                    if (min_resistant > max_susceptible) {
                        print min_resistant_row >> "'"$temp_file"'";
                        print max_susceptible_row >> "'"$temp_file"'";
                    }
                }
            
            #  If there is only data for one of the phenotypes, keep the highest value for 
            # Susceptible and the lowest for Resistant.
            } else if (has_resistant && !has_susceptible) {
            	if (min_resistant_row != "") {
            		print min_resistant_row >> "'"$temp_file"'"
            	} else if ( min_resistant_row == "") {
            		print min_resistant_row_2 >> "'"$temp_file"'"
            	}
            	
            } else if (has_susceptible && !has_resistant) {
            	if (max_susceptible_row != "") {
            		print max_susceptible_row >> "'"$temp_file"'"
            	} else if ( max_susceptible_row == "") {
            		print max_susceptible_row_2 >> "'"$temp_file"'"
            	}
            }
                  
        }
    }
    
}' "$input_file"


mv "$temp_file" "$output_file"
