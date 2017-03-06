#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [python, /home/bay001/projects/codebase/rbp-maps/maps/plot_density.py]

inputs:

  ip_bam:
    type: File
    format: http://edamontology.org/format_2572
    inputBinding:
      position: 1
      prefix: --ip
    label: "ip bam"
    doc: "ip bam"

  input_bam:
    type: File
    format: http://edamontology.org/format_2572
    inputBinding:
      position: 2
      prefix: --input
    label: "input bam"
    doc: "input bam"

  annotation:
    type: File
    format: http://edamontology.org/???
    inputBinding:
      position: 3
      prefix: -a
    label: "annotation"
    doc: "annotation bed rmats or miso file"

  annotation_type:
    type: string
    format: http://edamontology.org/???
    inputBinding:
      position: 4
      prefix: -t
    label: "annotation type"
    doc: "specify one of: bed rmats or miso"

  event:
    type: string
    format: http://edamontology.org/format_3320
    inputBinding:
      position: 5
      prefix: -e
    label: "event type"
    doc: "specify the type of splicing map to make (default se)"

arguments: [
  "-o",
  $(inputs.ip_bam.nameroot).svg
  ]

outputs:

  output_svg:
    type: File
    format: http://edamontology.org/format_3604
    outputBinding:
      glob: $(inputs.ip_bam.nameroot).svg
    label: ""
    doc: "rbp map svg file"
    
  output_raw_ip_density_map:
    type: File
    outputBinding:
      glob: $(inputs.ip_bam.nameroot).