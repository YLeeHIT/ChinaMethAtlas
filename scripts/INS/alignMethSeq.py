import pandas as pd
import os
import subprocess
from collections import defaultdict

minDepth = 3
minCG = 5

def average_meth(positions_values):
    result = []                             # 创建一个结果列表，用于存储每个位置的均值
    position_dict = defaultdict(list)       # 使用 defaultdict 来存储每个位置的值列表
    global minDepth
    for pos, value in positions_values:
        position_dict[pos].append(value)
    for pos in sorted(position_dict.keys()):
        values = position_dict[pos]
        # 至少需要3个C才可以进行计算
        print(f"pos is {pos} and value len is {len(values)}")
        if len(values) >= minDepth:
            avg_value = sum(values) / len(values)  # 计算均值
            result.append((pos, avg_value, len(values)))
    return result

def average_sites(position_values):
    averMeth = 0
    values = 0
    count = 0
    for pos, value, num in position_values:
        values += value
        count += 1
    averMeth = values/count if count > 0 else 0
    return averMeth


def sum_result(inputFile, outFile):
    cmd = f"cat {inputFile} >> {outFile}"
    # 执行命令
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("summary results")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running summary: {e.stderr.decode()}")


def insStepTwo(fm,fa):
    if os.path.getsize(fm) and os.path.getsize(fa):
        fmfile = pd.read_csv(fm,header=None, sep='\t')
        fafile = pd.read_csv(fa,header=None, sep='\t')
        posMeth = []
        for i in range(int(len(fafile)/2)):
            faReadName = fafile.iloc[2*i,0]
            faRead = fafile.iloc[2*i+1,0]
            faLen = len(faRead)
            fmReadNanme = fmfile.iloc[4*i,0]
            fmRead = fmfile.iloc[4*i+1,0]
            fmMeth = fmfile.iloc[4*i+3,0]
            fmLen = len(fmMeth)
            if faReadName == fmReadNanme:
                print(f"fm file is {fmRead}")
                print(f"fm name is {fmReadNanme}")
                faIndex = 0             # 带"_"间隔的 abpoa r1 得到文件
                fmIndex = 0             # 带meth level 的abpoa r0 得到文件，类似fq格式
                for i in range(faLen):
                    if faRead[i] != "-":
                        fmCG = fmRead[fmIndex:fmIndex+2]
                        if fmIndex + 1 < fmLen and fmCG == "CG":
                            posCG = faIndex + 1             # 输出文件坐标从1开始，不是0，方便substr提取
                            faCG = faRead[faIndex:faIndex+2]
                            methCG = int(fmMeth[fmIndex])
                            print(f"{posCG} and {fmCG} and {faCG} and {methCG}")
                            posMeth.append((posCG,methCG))
                        fmIndex += 1
                    faIndex += 1
        methResult = average_meth(posMeth)
        print(methResult)

        # 构造输出文件
        fileName = os.path.basename(fm).replace(".fm","").split("_")
        sampleName = fileName[0]
        chromosome = fileName[1]
        position = int(fileName[2])
        methLevel = average_sites(methResult)

        outDir = os.path.join(os.path.dirname(fm),"meth")
        os.makedirs(outDir, exist_ok=True)  # 创建输出目录
        outFile = os.path.join(outDir,os.path.basename(fm).replace(".fm",".ins"))
        print(f"Writting to {outFile}")
        with open(outFile, "w") as f:
            for pos, value, num in methResult:
                f.write(f"{pos}\t")   
                f.write(f"{value}\t")   
                f.write(f"{num}\n")   
        f.close()

        outDirResult = os.path.join(os.path.dirname(fm),"result")
        os.makedirs(outDirResult, exist_ok=True)  # 创建输出目录
        tmpName = ".".join([sampleName,"ins"])
        outFile2 = os.path.join(outDirResult,tmpName)
        print(f"Writting to {outFile2}")

        with open(outFile2, "w") as f:
            if len(methResult) >= minCG:
                f.write(f"{chromosome}\t")
                f.write(f"{position}\t")
                f.write(f"{position+1}\t")
                f.write(f"{methLevel:.4f}\t")
                f.write(f"{sampleName}\n")
        f.close()
        # 写到汇总文件中，记得去重
        outSummary = '/'.join(fm.split("/")[:-2])
        outResult = os.path.join(outSummary,".".join([sampleName,"meth"]))
        sum_result(outFile2,outResult)
    else:
        print(f"{fm} is empty") 
if __name__ == '__main__':
    fm = "/home/user/liyang/shell/DNAmeth/python/newSV/INS/data/xizang/100_DNA21SC026816-B/test/out/DNA21SC026816-B_chr1_710579.fm"
    fa = "/home/user/liyang/shell/DNAmeth/python/newSV/INS/data/xizang/100_DNA21SC026816-B/test/out/DNA21SC026816-B_chr1_710579_fa.fa"
    insStepTwo(fm,fa)
