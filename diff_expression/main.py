import kProcessor as kp


def differntialExpression(genes_file, samplesInput, controlInput, outputFilename):
    nSamples = len(samplesInput)
    nControl = len(controlInput)
    allDatasets = nSamples + nControl

    kFrames = list()  # kp.kFramesVector
    allSamples = list()
    allSamples = samplesInput + controlInput

    requiredIndices = list()

    index = 0

    count_col_name = "count"

    for filename in allSamples:
        currentFrame = kp.kDataFrameMQF()
        kp.loadFromKMC(currentFrame, filename)
        kp.createCountColumn(currentFrame)
        print(f"Loading {filename} kmers: {currentFrame.size()}")
        totalCount = kp.aggregate_count(currentFrame, count_col_name)
        currentFrame = kp.transform_normalize(currentFrame, count_col_name, totalCount)
        kFrames.append(currentFrame)

    kSize = int()
    if len(kFrames):
        kSize = kFrames[0].getkSize()

    chunkSize = 1000
    genesFrame = kp.kDataFrameMQF(kSize)
    kp.index(genesFrame, {"kSize": kSize}, genes_file, chunkSize, f"{genes_file}.names")
    kp.createColorColumn(genesFrame)
    kFrames.append(genesFrame)
    requiredIndices.append(len(kFrames)-1)
    colorColumn = f"color.{len(kFrames)-1}"
    print(f"Load {genes_file} kmers: {kFrames[-1].size()}")
    res = kp.innerJoin(kFrames, requiredIndices)
    res = kp.filter_zeroCounts(res, allDatasets, count_col_name)
    foldChange_col_name = "foldChange"
    kp.transform_foldchange(res, foldChange_col_name, count_col_name, nSamples, allDatasets)
    foldChangeByGene = kp.aggregate_foldChangeByGene(res, colorColumn)

    print(foldChangeByGene)


# Inputs
input_genes_file = str()
input_KMC_sample_files = list()
input_KMC_control_files = list()
output_filename = str()

differntialExpression(input_genes_file, input_KMC_sample_files, input_KMC_control_files, output_filename)
