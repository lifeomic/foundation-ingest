cwlVersion: v1.0
class: Workflow
hints:
  ResourceRequirement:
    coresMin: 1
    coresMax: 1
    ramMin: 3GB
    ramMax: 3GB
inputs:
  xmlFile: File
  reportFile: File
  source: string
  reference: string

outputs:
  copynumber:
    type: File
    outputSource: process_xml/copynumber
  pdf:
    type: File
    outputSource: process_xml/pdf
  structural:
    type: File
    outputSource: process_xml/structural
  tmp:
    type: File
    outputSource: process_xml/tmp
  vcf:
    type: File
    outputSource: normalize_vcf/normalized_vcf
  yml:
    type: File
    outputSource: generate_yml/yml

steps:
  process_xml:
    in:
      xmlFile: xmlFile
      reportFile: reportFile
      source: source
      LIFEOMIC_GENOME_REFERENCE: reference
    out:
      [copynumber, vcf, pdf, structural, tmp]
    run:
      class: CommandLineTool
      hints:
        DockerRequirement:
          dockerPull: lifeomic/foundation-ingest:2.3.0
      inputs:
        xmlFile:
          type: File
          inputBinding:
            prefix: -x
            position: 1
        reportFile:
          type: File
          inputBinding:
            prefix: -t
            position: 2
        source:
          type: string
          inputBinding:
            prefix: -s
            position: 3
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
            glob: '/tmp/**/*.report.pdf'
        structural:
          type: File
          outputBinding:
            glob: '/tmp/**/*.structural.csv'
        tmp:
          type: File
          outputBinding:
            glob: '/tmp/**/*.ga4gh.tmp'
  normalize_vcf:
    in:
      vcf:
        source: process_xml/vcf
      LIFEOMIC_GENOME_REFERENCE: reference
    out:
      [normalized_vcf]
    run:
      class: CommandLineTool
      hints:
        DockerRequirement:
          dockerPull: lifeomic/kopis-task-vtools:1.1.0
      arguments: ["vt-combo", "-r", "/tmp/reference/GRCh37.fa.gz"]
      inputs:
        vcf:
          type: File
          inputBinding:
            prefix: -i
            position: 1
        LIFEOMIC_GENOME_REFERENCE:
          type: string
      outputs:
        normalized_vcf:
          type: File
          outputBinding:
            glob: '/tmp/**/*.nrm.vcf'
  generate_yml:
    in:
      nrm_vcf:
        source: normalize_vcf/normalized_vcf
      yml_tmp:
        source: process_xml/tmp
    out:
      [yml]
    run:
      class: CommandLineTool
      hints:
        DockerRequirement:
          dockerPull: lifeomic/tmp-to-yml-ingest:1.0.0
      inputs:
        nrm_vcf:
          type: File
        yml_tmp:
          type: File
      outputs:
         yml:
          type: File
          outputBinding:
            glob: '/tmp/**/*.ga4gh.yml'

s:license: https://opensource.org/licenses/MIT
s:author:
  - class: s:Organization
    s:email: phc-clinical-research@lifeomic.com
    s:name: LifeOmic Clinical Research Team

$namespaces:
 s: https://schema.org/

$schemas:
 - https://schema.org/version/latest/schemaorg-current-http.rdf

doc: |
  Workflow that converts a Foundation XML test file into omics files to be
  ingested into the LifeOmic Precision Health Cloud (https://lifeomic.com/products/precision-health-cloud/).