import yaml
import glob
import re

chrs = [
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    "X"
]

with open("pog_libs.yaml") as f:
    pog_libs = yaml.load(f, Loader = yaml.FullLoader)

#print(pog_libs)

reverso = {}
for case in pog_libs.keys():
    reverso[case] = {}
    #print(pog_libs[case])
    for sample in pog_libs[case].keys():
#        print(sample)
         glob_res = glob.glob(pog_libs[case][sample]["merge"] + "/*dupsFlagged.bam")
         if not glob_res:
            continue
         reverso[case][pog_libs[case][sample]["name"]] = glob_res[0]

## Build input list

in_list = []

id_dict = pog_libs

cases=id_dict.keys()

blacklist=["A79702","A79703","A79704"]

#print(cases)

comps={}
bams={}
for case in cases:
    #print(id_dict[case])
    if id_dict[case].get("Normal_1") and id_dict[case].get("Diseased_1"):
#        print(case)
        norm = id_dict[case]["Normal_1"]
        bams[norm["name"]] = glob.glob(norm["merge"] + "/*dupsFlagged.bam")
        if len(bams[norm["name"]]) == 0:
            continue
        lib_iter = 1
        bams[norm["name"]] = bams[norm["name"]][0]
        tums = []
#        print(id_dict[case])
#        print("Diseased_" + str(lib_iter))
#        print(id_dict[case].get("Diseased_" + str(lib_iter)))
        while id_dict[case].get("Diseased_" + str(lib_iter)):
#            print(id_dict[case]["Diseased_" + str(lib_iter)])
            bam_glob = glob.glob(id_dict[case]["Diseased_" + str(lib_iter)]["merge"] + "/*dupsFlagged.bam")
            if len(bam_glob) == 0:
                lib_iter += 1
                continue
            if os.path.islink(bam_glob[0]) and not os.path.exists(bam_glob[0]):
                lib_iter += 1
                continue
            if id_dict[case]["Diseased_" + str(lib_iter)]["name"] in blacklist:
                lib_iter += 1
                continue
            bams[id_dict[case]["Diseased_" + str(lib_iter)]["name"]] = bam_glob[0]
            tums.append(id_dict[case]["Diseased_" + str(lib_iter)]["name"] + "_" + norm["name"])
            lib_iter += 1
        comps[case] = tums
        #print(tums)
        #libs = list(id_dict[case].keys())
        #libs = re.
        #print(libs)
in_list = []
for comp in comps.keys():
    in_list.extend(expand("data/{case}/{comp}/cna.txt", case = comp, comp = comps[comp]))


counter = 0
chunk = 1
chunk_dict = {}
chunk_list = []
for item in in_list:
    if counter == 49:
        counter = 0
        chunk_dict["chunk" + str(chunk)] = chunk_list
        chunk = chunk + 1
        chunk_list = []
    chunk_list.append("chunks/chunk" + str(chunk) + "/" + item)
    counter = counter + 1


in_list = chunk_dict.values()

chunks = chunk_dict.keys()
chunks = [int(re.sub("chunk([0-9]+)", "\\1", ch)) for ch in chunks]
chunk = "chunks/chunk" + str(max(chunks))
chunk_output = chunk + "/chunk_complete.txt"



#chunk_outputs = ["chunk" + str(num) for num in range(1, chunk)]
#chunk_outputs = "_".join(chunk_outputs)
#chunk_outputs = "chunk0_" + chunk_outputs 



print(chunk_output)

#completeness_dict = {}
#for chunk in chunk_dict.keys():
#    chunk_n = int(re.sub("chunk([0-9]+)", "\\1", chunk))
#    if chunk_n > 0:
#        chunk_n = ["chunk" + str(num) for num in range(chunk_n)]
#        chunk_n = "_".join(chunk_n)
#    elif chunk_n == 0:
#        chunk_n = "chunk0"
#    cases = chunk_dict[chunk]
#    for case in cases:
#        comp = re.sub(".*/([^/]*_[^/]*)/touched.txt", "\\1", case)
#        libs = comp.split("_")
#        for lib in libs:
#            completeness_dict[lib] = chunk_n

#print(completeness_dict)

def chunk_function(chunk):
    chunk = int(re.sub("chunk([0-9]+)", "\\1", chunk))
    chunk = "chunks/chunk" + str(chunk - 1)
    return(chunk + "/chunk_complete.txt")

localrules: all, norm_config, mk_tmp, gl_rm, init_chunk, run_chunk, seg_config

rule all:
    input:
        chunk_output

wildcard_constraints:
    case = "[^/]*",
    lib = "[^/]*",
    chunk = "[^/]*",
    lib1 = "[^/]*",
    lib2 = "[^/]*"

rule seq:
    input:
        lambda wildcards: bams[wildcards.lib],
        lambda wildcards: chunk_function(wildcards.chunk)
    output:
        temp(expand("chunks/{{chunk}}/data/{{case}}/{{lib}}/{chrs}.seq", chrs = chrs))
    resources: cpus=1, mem_mb=7900
    shell:
        "samtools/samtools view -U BWA,chunks/{wildcards.chunk}/data/{wildcards.case}/{wildcards.lib}/,N,N -q 20 {input[0]}"

rule norm_config:
    output:
        "chunks/{chunk}/data/{case}/{lib}/NORM_CONFIG.cfg"
    resources: cpus=1, mem_mb=7900
    shell:
        "python norm_config.py {wildcards.lib} chunks/{wildcards.chunk}/data/{wildcards.case}/{wildcards.lib}/ {output}"

rule mk_tmp:
    output:
        temp(directory("chunks/{chunk}/tmp_dir/{case}/{lib}/"))
    resources: cpus=1, mem_mb=7900
    shell:
        "mkdir -p {output}"

rule gl_rm:
    input:
        expand("chunks/{{chunk}}/data/{{case}}/{{lib}}/{chrs}.seq", chrs = chrs)
    output:
        "chunks/{chunk}/data/{case}/{lib}/gl_rm.status"
    resources: cpus=1, mem_mb=7900
    shell:
        "rm chunks/{wildcards.chunk}/data/{wildcards.case}/{wildcards.lib}/GL*.seq; touch {output}"

rule norm:
    input:
        expand("chunks/{{chunk}}/data/{{case}}/{{lib}}/{chrs}.seq", chrs = chrs),
        "chunks/{chunk}/data/{case}/{lib}/NORM_CONFIG.cfg",
        "chunks/{chunk}/tmp_dir/{case}/{lib}/",
        "chunks/{chunk}/data/{case}/{lib}/gl_rm.status"
    output:
        expand("chunks/{{chunk}}/data/{{case}}/{{lib}}/{chrs}.norm.bin", chrs = chrs),
        temp("chunks/{chunk}/data/{case}/{lib}/tmpfile")
    resources: cpus=1, mem_mb=7900
    shell:
        "./bicseq_install/NBICseq-norm_v0.2.4/NBICseq-norm.pl chunks/{wildcards.chunk}/data/{wildcards.case}/{wildcards.lib}/NORM_CONFIG.cfg --tmp=chunks/{wildcards.chunk}/tmp_dir/{wildcards.case}/{wildcards.lib}/ -l 150 chunks/{wildcards.chunk}/data/{wildcards.case}/{wildcards.lib}/tmpfile"

rule init_chunk:
    output:
        "chunks/chunk0/chunk_complete.txt"
    resources: cpus=1, mem_mb=7900
    shell:
        "touch {output}"

rule seg_config:
    output:
        "chunks/{chunk}/data/{case}/{lib1}_{lib2}/SEG_CONFIG.cfg"
    resources: cpus=1, mem_mb=7900
    shell:
        "python seg_config.py chunks/{wildcards.chunk}/data/{wildcards.case}/{wildcards.lib2}/ chunks/{wildcards.chunk}/data/{wildcards.case}/{wildcards.lib1}/ {output}"

rule cna:
    input:
        ancient(expand("chunks/{{chunk}}/data/{{case}}/{{lib1}}/{chrs}.norm.bin", chrs = chrs)),
        ancient(expand("chunks/{{chunk}}/data/{{case}}/{{lib2}}/{chrs}.norm.bin", chrs = chrs)),
        "chunks/{chunk}/data/{case}/{lib1}_{lib2}/SEG_CONFIG.cfg"
    output:
        "chunks/{chunk}/data/{case}/{lib1}_{lib2}/cna.txt"
    resources: cpus=1, mem_mb=7900
    shell:
        "./bicseq_install/NBICseq-seg_v0.7.2/NBICseq-seg.pl chunks/{wildcards.chunk}/data/{wildcards.case}/{wildcards.lib1}_{wildcards.lib2}/SEG_CONFIG.cfg {output}"

rule run_chunk:
    input:
        lambda wildcards: chunk_function(wildcards.chunk),
        lambda wildcards: chunk_dict[wildcards.chunk]
    output:
        "chunks/{chunk}/chunk_complete.txt"
    resources: cpus=1, mem_mb=7900
    shell:
        "touch {output}"
