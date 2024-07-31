cd /fh/fast/gilbert_p/fg_data/flowplyr
ml fhR/4.2.0-foss-2021b
Rscript startup.R \
--run_name 'unit_test' \
--base_dir '/fh/fast/gilbert_p/fg_data/ics_test/flowplyr_test_data/' \
--fcm08_path '/fh/fast/gilbert_p/fg_data/ics_test/flowplyr_test_data/fcm08_simulated_format.txt' \
--fcs_folder_path '/fh/fast/gilbert_p/fg_data/ics_test/flowplyr_test_data/3639-L-ICS PROFICIENCY 3/' \
--xml_file_path '/fh/fast/gilbert_p/fg_data/ics_test/flowplyr_test_data/3639-L-ICS PROFICIENCY 3/3639-L-ICS PROFICIENCY 3 FJ.xml' \
--parent_gate "/Time/S/14-/Lv/L/3+/4+" \
--stim_exclusion_terms 'phactrl,sebctrl,posctrl' \
--functional_markers '154+,GzB+,IFNg+,IL2+,IL4+,IL17a+' \
--xml_keywords '$FIL,Stim,Sample Order,EXPERIMENT NAME,Replicate'



readr::read_tsv('/fh/fast/gilbert_p/fg_data/ics_test/flowplyr_test_data/fcm08_simulated.txt')


