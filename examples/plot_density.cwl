#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [plot_density.py]

inputs:

  ip_bam:
    type: File
    format: http://edamontology.org/format_2572
    inputBinding:
      position: 1
      prefix: --ipbam
    label: "ip bam"
    doc: "ip bam"

  input_bam:
    type: File
    format: http://edamontology.org/format_2572
    inputBinding:
      position: 2
      prefix: --inputbam
    label: "input bam"
    doc: "input bam"

  annotations:
    type: File[]
    format: http://edamontology.org/???
    inputBinding:
      position: 3
      prefix: --annotations
    label: "annotations"
    doc: "annotation bed rmats or miso file(s)"

  annotation_type:
    type: string[]
    format: http://edamontology.org/???
    inputBinding:
      position: 4
      prefix: --annotation_type
    label: "annotation type"
    doc: "specify one of: bed rmats or miso for each annotation"

  event:
    type: string
    format: http://edamontology.org/format_3320
    inputBinding:
      position: 5
      prefix: --event
    label: "event type"
    doc: "specify the type of splicing map to make (default se)"

  exon_offset:
    type: int
    format: http://edamontology.org/format_3320
    inputBinding:
      position: 6
      prefix: --exon_offset
    label: "exon offset"
    doc: "specify exon offset"

  intron_offset:
    type: int
    format: http://edamontology.org/format_3320
    inputBinding:
      position: 7
      prefix: --intron_offset
    label: "intron offset"
    doc: "specify intron offset"

  scale:
    type: boolean
    default: False
    inputBinding:
      position: 8
      prefix: --scale
    label: "scaled"
    doc: "specify if events should be scaled"

  chrom_sizes:
    type: File
    format: http://edamontology.org/???
    inputBinding:
      position: 9
      prefix: --chrom_sizes
    label: "chrom sizes"
    doc: "chrom sizes file from ucsc"

  unflipped:
    type: boolean
    default: False
    inputBinding:
      position: 10
      prefix: --unflipped
    label: "unflipped"
    doc: "specify if the *pos.bw and *neg.bw are actually flipped strands"

arguments: [
  "--output",
  $(inputs.ip_bam.nameroot).$(inputs.event).svg
  ]

outputs:

  output_svg:
    type: File
    format: http://edamontology.org/format_3604
    outputBinding:
      glob: $(inputs.ip_bam.nameroot).$(inputs.event).svg
    label: ""
    doc: "rbp map svg file"

  output_all_means_map:
    type: File[]
    outputBinding:
      glob: $(inputs.ip_bam.nameroot).*.txt

  output_all_density_map:
    type: File[]
    outputBinding:
      glob: $(inputs.ip_bam.nameroot).*.csv