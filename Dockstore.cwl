cwlVersion: v1.0
class: CommandLineTool
id: foundation-ingest
label: foundation-ingest
hints:
  DockerRequirement:
    dockerPull: lifeomic/foundation-ingest:1.0.2

inputs:
  xmlFile:
    type: File
    inputBinding:
      prefix: -x
      position: 1
  source:
    type: string
    inputBinding:
      prefix: -s
      position: 2
  LIFEOMIC_GENOME_REFERENCE:
    type: string

outputs:
  copynumber:
    type: File
    outputBinding:
      glob: '/tmp/**/*.copynumber.csv'
  vcf:
    type: File
    outputBinding:
      glob: '/tmp/**/*.vcf'
  pdf:
    type: File
    outputBinding:
      glob: '/tmp/**/*.pdf'
  structural:
    type: File
    outputBinding:
      glob: '/tmp/**/*.structural.csv'
  tmp:
    type: File
    outputBinding:
      glob: '/tmp/**/*.ga4gh.tmp'

doc: |
  This tool converts a Foundation XML file into omic files that can be
  ingested into the LifeOmic Precision Health Cloud (https://lifeomic.com/products/precision-health-cloud/).