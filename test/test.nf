params.data_root = false

workflow {
    root = file(params.data_root)
    t1_channel = Channel.fromFilePairs("$root/**/*_t1.nii.gz", size: 1, flat: true)
    //t1_channel.view()
    seg_channel = Channel.fromFilePairs("$root/**/*{wm,gm,csf}_mask.nii.gz", size: 3, flat: true)
    seg_channel = t1_channel.join(seg_channel, remainder: true).map{ it.size() > 3 ? it.subList(0, it.size() - 1) : [it[0], "", "", ""] }
    seg_channel.view()
    seg_channel.count().view()
}