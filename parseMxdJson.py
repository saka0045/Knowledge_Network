import json
import argparse
import os
from bs4 import BeautifulSoup


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inputJson", dest="input_path", required=True,
        help="Full path to the input MXD json file"
    )
    parser.add_argument(
        "-t", "--templateFile", dest="template_path", required=True,
        help="Full path to the template Knowledge Network upload csv"
    )
    parser.add_argument(
        "-o", "--outFile", dest="output_path", required=True,
        help="Full path to the output csv file"
    )
    parser.add_argument(
        "-b", "--bedFile", dest="bed_path", required=True,
        help="Path to the hg38 bed file"
    )
    parser.add_argument(
        "-f", "--refBaseFile", dest="ref_base_path", required=True,
        help="Path to the file with reference base in hg38"
    )

    args = parser.parse_args()

    input_path = os.path.abspath(args.input_path)
    template_path = os.path.abspath(args.template_path)
    output_path = os.path.abspath(args.output_path)
    bed_path = os.path.abspath(args.bed_path)
    ref_base_path = os.path.abspath(args.ref_base_path)

    open_json = open(input_path, "r")
    template_file = open(template_path, "r")
    upload_file = open(output_path, "w")
    bed_file = open(bed_path, "r")
    ref_base_file = open(ref_base_path, "r")

    # Open json file
    json_file = json.load(open_json)

    # Go to the "results" section of the json file
    json_data = json_file["results"]

    # Copy the header from the template file
    for line in template_file:
        line = line.rstrip()
        upload_file.write(line)
        upload_file.write("\n")

    # Parse the json file
    for entry in range(0, len(json_data)):
        mxd_data = extract_mxd_data(json_data, entry)

        # Write the above information into the upload file
        # Skip chr X, it failed the hg38 lift over for now
        if mxd_data[12] == "x":
            continue
        else:
            # Replace the hg19 coordinates with hg38
            hg38_coordinates = bed_file.readline()
            hg38_coordinates = hg38_coordinates.rstrip()
            hg38_coordinates_list = hg38_coordinates.split("\t")
            hg38_start = str(int(hg38_coordinates_list[1]) + 1)
            hg38_stop = str(int(hg38_coordinates_list[2]) + 1)
            # Get the ref base for hg38
            ref_base_line = ref_base_file.readline()
            ref_base_line = ref_base_line.rstrip()
            ref_base_items = ref_base_line.split("\t")
            hg38_ref_base = ref_base_items[1]
            # If kn_variant_type is Variant, only use start position
            if mxd_data[7] == "Variant":
                mxd_data[13] = hg38_start
                hg19_ref_base = mxd_data[14]
                hg19_alt_base = mxd_data[15]
                # Check to see if the ref bases are the same across the two genomes
                if hg38_ref_base == hg19_ref_base:
                    upload_file.write(",".join(mxd_data))
                    upload_file.write("\n")
                # If the hg19 alt base became the hg38 ref base, skip the variant
                elif hg38_ref_base == hg19_alt_base:
                    print("Variant " + mxd_data[0] + " hg19 alt base: " + hg19_alt_base + " became the hg38 ref base: "
                          + hg38_ref_base)
                    print("Skipping variant")
                    continue
                elif hg38_ref_base != hg19_ref_base:
                    # if the ref base is different the lengths don't match, ignore for now
                    if len(hg38_ref_base) != hg19_ref_base:
                        upload_file.write(",".join(mxd_data))
                        upload_file.write("\n")
                    # if the ref base is different and the lengths do match, replace the hg19 ref with hg38 ref
                    elif len(hg38_ref_base) == len(hg19_ref_base):
                        print("Variant " + mxd_data[0] + " hg19 ref base " + hg19_ref_base +
                              " is different from hg38 ref base: " + hg38_ref_base)
                        mxd_data[14] = hg38_ref_base
                        upload_file.write(",".join(mxd_data))
                        upload_file.write("\n")
            # If kn_variant_type is SV, use start and stop
            elif mxd_data[7] == "SV":
                mxd_data[13] = (hg38_start + ":" + hg38_stop)
                upload_file.write(",".join(mxd_data))
                upload_file.write("\n")

    open_json.close()
    upload_file.close()
    template_file.close()
    bed_file.close()
    ref_base_file.close()


def extract_mxd_data(json_data, entry):
    """
    Extracts the MXD data and converts it to Knowledge Network uploadable format.
    Used within the loop to iterate through all queries for the json file
    :param json_data:
    :param entry:
    :return:
    """
    # Create an empty list to store the below variables
    mxd_data_for_entry = []
    # Get the variant type
    try:
        coordinateType = json_data[entry]["transientVariantDetails"]["coordinates"]["hg19"]["coordinateType"]
        # MXD seems to be inconsistent on naming this variantSubtype or some don't have any
        if coordinateType == "simple":
            mxd_variant_type = "snv"
        elif coordinateType == "cnv":
            mxd_variant_type = "cnv"
    except TypeError:
        coorfinateType = "snv"
    # print(coordinateType + " and " + mxd_variant_type)

    referenceId = json_data[entry]["id"]
    # print(referenceId)
    pathogenicity = json_data[entry]["conclusion"]
    # Capitalize each word
    pathogenicity = pathogenicity.title()
    inheritanceMode = "Not Specified"
    confidence = ""
    penetrance = ""
    mechanism = ""
    try:
        html_text = json_data[entry]["note"]
        # Convert the html to a text
        associationCuratorSummary = BeautifulSoup(html_text, features="lxml")
        associationCuratorSummary = associationCuratorSummary.get_text()
        # Escape any double quotes in the comment
        associationCuratorSummary = associationCuratorSummary.replace("\"", "\"\"")
        # Escape any "," in the comment by surrounding the entire comment with double quotes
        associationCuratorSummary = "\"" + associationCuratorSummary + "\""
        # Replace new line characters
        associationCuratorSummary = associationCuratorSummary.replace("\n", "")
    except KeyError:
        associationCuratorSummary = ""
    except TypeError:
        associationCuratorSummary = ""
    # Type and curation data depends on the variant_type
    if mxd_variant_type == "snv":
        kn_variant_type = "Variant"
        curationLevel = "Nucleotide"
    elif mxd_variant_type == "cnv":
        kn_variant_type = "SV"
        curationLevel = "Position"
    mutationClass = "Not Specified"
    multipleBiomarker = ""
    vid = ""
    chromosome = json_data[entry]["transientVariantDetails"]["coordinates"]["hg19"]["reference"]
    # SNVs will only have the start position and also has ref and alt
    if mxd_variant_type == "snv":
        position = str(json_data[entry]["transientVariantDetails"]["coordinates"]["hg19"]["start"])
        ref = json_data[entry]["transientVariantDetails"]["coordinates"]["hg19"]["referenceBases"]
        alt = json_data[entry]["transientVariantDetails"]["coordinates"]["hg19"]["alternateBases"]
    # For CNVs, it needs start:stop and the ref and alt will be empty
    elif mxd_variant_type == "cnv":
        start_position = str(json_data[entry]["transientVariantDetails"]["coordinates"]["hg19"]["start"]["position"])
        stop_position = str(json_data[entry]["transientVariantDetails"]["coordinates"]["hg19"]["end"]["position"])
        position = (start_position + ":" + stop_position)
        ref = ""
        alt = ""
    aminoAcid = ""
    codon = ""
    exonNumber = ""
    # This is a list of genes
    geneSymbol = json_data[entry]["gene"]
    # Capitalize the genes
    geneSymbol = [gene.replace(gene, gene.upper()) for gene in geneSymbol]
    transcriptId = json_data[entry]["additionalDetails"]["transcript"]
    multipleDiseases = ""
    # Using the generic "Disease (disorder)" for disease name since there are no entries in MXD
    diseaseName = "Disease (disorder)"
    diseaseOntologyGroup = "SNOMEDCT"
    diseaseOntologyGroupId = "64572001"
    evidenceType = "Miscellaneous"
    supportType = "Not Selected"
    evidenceStrength = "Not Selected"
    criteria = ""
    evidenceCuratorSummary = ""
    evidenceCategory = ""
    source = ""
    content = ""
    pubmedId = ""
    url = ""
    organization = ""
    licenseType = ""
    openAccess = ""
    mxd_data_for_entry.extend((referenceId, pathogenicity, inheritanceMode, confidence, penetrance, mechanism,
                     associationCuratorSummary, kn_variant_type, curationLevel, mutationClass, multipleBiomarker,
                     vid, chromosome, position, ref, alt, aminoAcid, codon, exonNumber, " ".join(geneSymbol),
                     transcriptId, multipleDiseases, diseaseName, diseaseOntologyGroup, diseaseOntologyGroupId,
                     evidenceType, supportType, evidenceStrength, criteria, evidenceCuratorSummary, evidenceCategory,
                     source, content, pubmedId, url, organization, licenseType, openAccess))

    return mxd_data_for_entry


if __name__ == "__main__":
    main()
