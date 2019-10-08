import json

open_json = open("/Users/m006703/Illumina/Knowledge_Network/files/nmpan.json", "r")
upload_file = open("/Users/m006703/Illumina/Knowledge_Network/files/KN_Upload3.csv", "w")
template_file = open("/Users/m006703/Illumina/Knowledge_Network/files/KN_Upload_Template.csv", "r")
'''
# This portion was to deal with the first set of MXD query from Eric

# Make the MXD download file a valid json file
json_without_comments.write("[\n")
for line in open_json:
    if line.startswith("/*"):
        continue
    else:
        line = line.replace("ObjectId(", "")
        line = line.replace(")", "")
        line = line.replace("ISODate(", "")
        line = line.replace("NumberLong(", "")
        json_without_comments.write(line)

json_without_comments.write("\n]")
json_without_comments.close()

json_without_comments = open("/Users/m006703/Illumina/Knowledge_Network/Modified_Curations.json", "r")

json_data = json.load(json_without_comments)
'''
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
    # Create an empty list to store the below variables
    mxd_data = []

    referenceId = json_data[entry]["id"]
    print(referenceId)
    pathogenicity = json_data[entry]["conclusion"]
    # Capitalize each word
    pathogenicity = pathogenicity.title()
    inheritanceMode = "Not Specified"
    confidence = ""
    penetrance = ""
    mechanism = ""
    try:
        associationCuratorSummary = json_data[entry]["note"]
        # Can't pull comments with the latest example from Eric, disabling this for now
        associationCuratorSummary = ""
        # Escape any "," in the note
        #associationCuratorSummary = "\"" + associationCuratorSummary + "\""
    except KeyError:
        associationCuratorSummary = ""
    except TypeError:
        associationCuratorSummary = ""
    # If MXD data has something other than SNV, the following two variables needs to be changed
    type = "Variant"
    curationLevel = "Nucleotide"
    mutationClass = "Not Specified"
    multipleBiomarker = ""
    vid = ""
    chromosome = json_data[entry]["transientVariantDetails"]["coordinates"]["hg19"]["reference"]
    # If the variant is SV, it will also need the stop position (start:stop)
    position = str(json_data[entry]["transientVariantDetails"]["coordinates"]["hg19"]["start"])
    ref = json_data[entry]["transientVariantDetails"]["coordinates"]["hg19"]["referenceBases"]
    alt = json_data[entry]["transientVariantDetails"]["coordinates"]["hg19"]["alternateBases"]
    aminoAcid = ""
    codon = ""
    exonNumber = ""
    # This is a list of genes
    geneSymbol = json_data[entry]["gene"]
    # Capitalize the genes
    geneSymbol = [gene.replace(gene, gene.upper()) for gene in geneSymbol]
    transcriptId = json_data[entry]["additionalDetails"]["transcript"]
    multipleDiseases = ""
    diseaseName = ""
    diseaseOntologyGroup = ""
    diseaseOntologyGroupId = ""
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

    mxd_data.extend((referenceId, pathogenicity, inheritanceMode, confidence, penetrance, mechanism,
                     associationCuratorSummary, type, curationLevel, mutationClass, multipleBiomarker,
                     vid, chromosome, position, ref, alt, aminoAcid, codon, exonNumber, " ".join(geneSymbol),
                     transcriptId, multipleDiseases, diseaseName, diseaseOntologyGroup, diseaseOntologyGroupId,
                     evidenceType, supportType, evidenceStrength, criteria, evidenceCuratorSummary, evidenceCategory,
                     source, content, pubmedId, url, organization, licenseType, openAccess))

    # Write the above information into the upload file
    upload_file.write(",".join(mxd_data))
    upload_file.write("\n")

open_json.close()
upload_file.close()
template_file.close()
