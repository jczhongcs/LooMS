#!/bin/bash
"""
This script is to highlight the flow of running WaterlooMS from RAW data acquisition to PSM output. 
"""



cd /home/jia/Documents/Code/waterlooms/dia_data_reading


mkdir data/PXD005573/

wget -P data/PXD005573/ ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2017/10/PXD005573/Fig1_MP-DIA-120min-120kMS1-10W30k-14dppp_MHRM_R01.raw
wget -P data/PXD005573/ ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2017/10/PXD005573/uniprot_sprot_2014-12-11_HUMAN_ISOFORMS.fasta


# Convert the data set, with centroiding
docker run -it -v $PWD/data/PXD005573:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:x64 wine msconvert --zlib --filter "peakPicking true 1-" --mzML /data/Fig1_MP-DIA-120min-120kMS1-10W30k-14dppp_MHRM_R01.raw

# Convert the file to be owned by the USER:GROUP of the host 
docker run -it -v $PWD/data/PXD005573:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:x64 chown $(id -u ${USER}):$(id -g ${USER}) Fig1_MP-DIA-120min-120kMS1-10W30k-14dppp_MHRM_R01.mzML



# TODO: Debug this with Johra
# Download the .jar release of the deBruijn decoy generator
wget -P data/PXD005573/ https://github.com/johramoosa/deBruijn/blob/master/deBruijn.jar?raw=true

mv data/PXD005573/deBruijn.jar?raw=true data/PXD005573/deBruijn.jar

cd data/PXD005573/
java -jar /home/jia/Documents/Code/waterlooms/dia_data_reading/data/PXD005573/deBruijn.jar --input /home/jia/Documents/Code/waterlooms/dia_data_reading/data/PXD005573/uniprot_sprot_2014-12-11_HUMAN_ISOFORMS.fasta

# Clean the temporary files from deBruijn
./clean_deBruijn.sh




# Copy the filename for the data over for a more weildy filename (optional)
mv data/PXD005573/Fig1_MP-DIA-120min-120kMS1-10W30k-14dppp_MHRM_R01.mzML data/PXD005573/dia_data.mzML

# Invoke WaterlooMS on the now converted mzMmL file
# Run the Maven lifecycle to generate a FatJar of WaterlooMSs
mvn clean install package






docker run -it -v $PWD:/data waterlooms:1.0 java -jar waterlooms.jar -detectionParams featuredetect.params -selectionParams featureselection.params









# MISC
cat ~/GH_TOKEN.txt | docker login docker.pkg.github.com -u JRWu --password-stdin


"""
For DEV you want to mount the entire directory and have access to the file(s)
However for PROD you want only the /src/main/python directory mounted





"""








