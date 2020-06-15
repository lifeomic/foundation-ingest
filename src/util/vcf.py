import re
import os

try:
    import pyhgvs as hgvs
    import pyhgvs.utils as hgvs_utils
    from pyfaidx import Fasta
except ImportError:
    Fasta = None


_COMP = dict(A="T", C="G", G="C", T="A", N="N", a="t", c="g", g="c", t="a", n="n")


def hgvs_2_vcf(
    variant_name, genes, functional_effect, cds_effect, position_value, strand, fasta
):
    if functional_effect in ["splice", "frameshift", "nonframeshift"]:
        return parse_splice(cds_effect, position_value, strand, fasta)
    else:
        try:
            return parse_hgvs(variant_name, fasta, genes)
        except:
            return parse_splice(cds_effect, position_value, strand, fasta)


def extract_vcf(variants, specimen_name, fasta, genes, output):
    with open("/tmp/unsorted.vcf", "w+") as vcf_file:
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write("##source=foundation-xml-fhir\n")
        vcf_file.write(f"##reference=file://{fasta}\n")
        vcf_file.write(
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
        )
        vcf_file.write(
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        )
        vcf_file.write(
            '##INFO=<ID=VENDSIG,Number=1,Type=String,Description="Vendor Significance">\n'
        )
        vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf_file.write(
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
        )
        vcf_file.write(
            '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Number of reads harboring allele (in order specified by GT)">\n'
        )
        vcf_file.write("##contig=<ID=chr1,length=248956422>\n")
        vcf_file.write("##contig=<ID=chr2,length=242193529>\n")
        vcf_file.write("##contig=<ID=chr3,length=198295559>\n")
        vcf_file.write("##contig=<ID=chr4,length=190214555>\n")
        vcf_file.write("##contig=<ID=chr5,length=181538259>\n")
        vcf_file.write("##contig=<ID=chr6,length=170805979>\n")
        vcf_file.write("##contig=<ID=chr7,length=159345973>\n")
        vcf_file.write("##contig=<ID=chr8,length=145138636>\n")
        vcf_file.write("##contig=<ID=chr9,length=138394717>\n")
        vcf_file.write("##contig=<ID=chr10,length=133797422>\n")
        vcf_file.write("##contig=<ID=chr11,length=135086622>\n")
        vcf_file.write("##contig=<ID=chr12,length=133275309>\n")
        vcf_file.write("##contig=<ID=chr13,length=114364328>\n")
        vcf_file.write("##contig=<ID=chr14,length=107043718>\n")
        vcf_file.write("##contig=<ID=chr15,length=101991189>\n")
        vcf_file.write("##contig=<ID=chr16,length=90338345>\n")
        vcf_file.write("##contig=<ID=chr17,length=83257441>\n")
        vcf_file.write("##contig=<ID=chr18,length=80373285>\n")
        vcf_file.write("##contig=<ID=chr19,length=58617616>\n")
        vcf_file.write("##contig=<ID=chr20,length=64444167>\n")
        vcf_file.write("##contig=<ID=chr21,length=46709983>\n")
        vcf_file.write("##contig=<ID=chr22,length=50818468>\n")
        vcf_file.write("##contig=<ID=chrX,length=156040895>\n")
        vcf_file.write("##contig=<ID=chrY,length=57227415>\n")
        vcf_file.write("##contig=<ID=chrM,length=16569>\n")
        vcf_file.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{specimen_name}\n"
        )

        status = {
            "known": "Pathogenic",
            "likley": "Likely_pathogenic",
            "unknown": "Uncertain_significance",
            "ambiguous": "other",
        }

        for variant_dict in variants:
            vendsig = status.get(variant_dict.get("@status", "unknown"))
            cds_effect = variant_dict["@cds-effect"].replace("&gt;", ">")
            transcript = variant_dict["@transcript"]
            functional_effect = variant_dict["@functional-effect"]
            strand = variant_dict["@strand"]
            position_value = variant_dict["@position"]
            dp = variant_dict["@depth"]
            af = variant_dict["@allele-fraction"]
            gt = "1/1" if float(variant_dict["@allele-fraction"]) > 0.9 else "0/1"
            alt = int(round(int(dp) * float(af)))
            ref = int(dp) - alt
            ad = f"{ref},{alt}"
            variant_name = f"{transcript}:c.{cds_effect}"

            chrom, offset, ref, alt = hgvs_2_vcf(
                variant_name,
                genes,
                functional_effect,
                cds_effect,
                position_value,
                strand,
                fasta,
            )
            vcf_file.write(
                f"{chrom}\t{offset}\t.\t{ref}\t{alt}\t.\tPASS\tDP={dp};AF={af};VENDSIG={vendsig}\tGT:DP:AD\t{gt}:{dp}:{ad}\n"
            )

    vcf_name = f'{output}/{specimen_name}/{specimen_name}.vcf'
    os.system(f'grep "^#" /tmp/unsorted.vcf > {vcf_name} && grep -v "^#" /tmp/unsorted.vcf | sort -V -k1,1 -k2,2n >> {vcf_name}')


def parse_hgvs(hgvs_name, fasta, genes):
    genome = Fasta(fasta, key_function=lambda x: f"chr{x}")

    with open(genes) as infile:
        transcripts = hgvs_utils.read_transcripts(infile)

    def get_transcript(name):
        return transcripts.get(name)

    return hgvs.parse_hgvs_name(hgvs_name, genome, get_transcript=get_transcript)


def getRevComp(seq):
    return "".join(_COMP[base] for base in reversed(seq))


def getSequence(genome, chrom, start, end):
    return str(genome[str(chrom)][start - 1 : end]).upper()


def parse_splice(cdsEffect, position, strand, fasta):
    genome = Fasta(fasta, key_function=lambda x: f"chr{x}")

    [chr, sPos] = position.split(":")
    startPos = int(sPos)
    if ">" in cdsEffect:
        mylist = cdsEffect.split(">")
        mylist[0] = re.sub(r"[A-Za-z]", "", mylist[0])
        mytype = "sub"
    elif ">" in cdsEffect:
        mylist = cdsEffect.split(">")
        mylist[0] = re.sub(r"[A-Za-z]", "", mylist[0])
        mytype = "sub"
    elif "delins" in cdsEffect:
        mylist = cdsEffect.split("delins")
        mytype = "delins"
    elif "del" in cdsEffect:
        mylist = cdsEffect.split("del")
        mytype = "del"
    elif "ins" in cdsEffect:
        mylist = cdsEffect.split("ins")
        mytype = "ins"
    elif "dup" in cdsEffect:
        mylist = cdsEffect.split("dup")
        mytype = "dup"
    else:
        raise ValueError("ERROR: not sure how to interpret [{}]")

    if strand == "-" and re.match(r"^[A-Za-z]+", mylist[1]):
        mylist[1] = getRevComp(mylist[1])
    myrange = mylist[0].split("_")
    if len(myrange) == 1:
        if mytype == "sub":
            ref = getSequence(genome, chr, startPos, startPos)
            return (chr, startPos, ref, mylist[1])
        elif mytype == "del":
            ref = getSequence(genome, chr, startPos, startPos + 1)
            return (chr, startPos, ref, ref[0])
        elif mytype == "dup":
            ref = getSequence(genome, chr, startPos, startPos + 1)
            return (chr, startPos, ref[0], ref)
        else:
            raise ValueError(
                f"ERROR: not a range value and not a substitution or deletion [{cdsEffect}]"
            )
    elif len(myrange) == 2:
        main = ["", ""]
        sub = ["", ""]
        for i in range(len(myrange)):
            if "+" in myrange[i]:
                [main[i], sub[i]] = re.split("\+", myrange[i])
            elif myrange[i].startswith("-"):
                main[i] = myrange[i]
                sub[i] = "0"
            elif "-" in myrange[i]:
                [main[i], sub[i]] = re.split("-", myrange[i])
                sub[i] = "-" + sub[i]
            else:
                main[i] = myrange[i]
                sub[i] = "0"

        if re.match(r"^[0-9]+", mylist[1]):
            seqlen = int(mylist[1])
        else:
            seqlen = len(mylist[1])

        if main[0].startswith("*") or main[1].startswith("*"):
            mylen = seqlen
        elif main[0] == main[1]:
            mylen = abs(int(sub[0]) - int(sub[1])) + 1
        else:
            if sub[0] != "0" and sub[1] != "0":
                if int(main[1]) - int(main[0]) == 1:
                    mylen = seqlen
                else:
                    raise ValueError(
                        "ERROR: main position has to match for both intronic else one of them must be an exon: "
                    )
            elif sub[0] == "0":
                mylen = abs(int(main[0]) - int(main[1])) + 1
                mylen += abs(int(sub[1]))
            elif sub[1] == "0":
                mylen = abs(int(main[0]) - int(main[1])) + 1
                mylen += abs(int(sub[0]))

        if mytype == "ins":
            if mylen != 2:
                raise ValueError(
                    f"ERROR: insertion but range is not 1 [{cdsEffect}]"
                )
            elif not re.match(r"^[A-Za-z]+", mylist[1]):
                ref = getSequence(genome, chr, startPos, startPos)
                mylist[1] = ref + ("N" * int(mylist[1]))
            else:
                ref = getSequence(genome, chr, startPos, startPos)
                mylist[1] = ref + mylist[1]
        elif mytype == "del":
            #    startPos = startPos - 1
            #    if mylen != seqlen:
            #        raise ValueError('ERROR: length of cds range does not match the given deleted sequence  [{}]'.format(cdsEffect))
            ref = getSequence(genome, chr, startPos, startPos + seqlen)
            mylist[1] = ref[0]
        elif mytype == "dup":
            if mylen != seqlen:
                raise ValueError(
                    f"ERROR: length of cds range does not match the given duplicated sequence  [{cdsEffect}]"
                )
            ref = getSequence(genome, chr, startPos, startPos + mylen)
            mylist[1] = ref
            ref = ref[0]
        elif main[0].startswith("*") or main[1].startswith("*"):
            raise ValueError(
                f"ERROR: insert+delete on 5UTR, unable to resolve sequence [{cdsEffect}]"
            )
        else:
            ref = getSequence(genome, chr, startPos, startPos + mylen - 1)

        return (chr, startPos, ref, mylist[1])
    else:
        raise ValueError("ERROR: syntax dont fit")

