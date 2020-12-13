#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def key_from_filename ( chan, split_char ) {
    return chan.map{
        def fname = it.getFileName().toString().split(split_char)
        [ fname[0..<(fname.size() > 2 ? 2 : 1)].join(split_char), it ]
    }
}

def group_channel_rep ( chan ) {
    return chan.groupTuple().map{
        [it[0]] + it.subList(1, it.size()).inject((1..it.size()).collect{ [] }) { sub, rep ->
            sub.eachWithIndex { its, i -> its.add(rep[i]) } ; return sub
        }
    }
}

def group_subject_reps ( dwi_channel, metadata_channel ) {
    return [group_channel_rep(dwi_channel), group_channel_rep(metadata_channel)]
}

def replace_dwi_file ( base_channel, dwi_channel ) {
    return dwi_channel.join(base_channel.map{ [it[0]] + it.subList(2, it.size()) })
}

OPT_FILE_VALUE = ""
OPT_CHANNEL = null

def opt_channel () {
    return OPT_CHANNEL
}

def join_optional ( base_channel, opt_channel ) {
    if ( opt_channel )
        return base_channel.join(opt_channel)
    else
        return base_channel.map{ it + [OPT_FILE_VALUE] }
}

def map_optional ( base_channel, opt_idx ) {
    return base_channel.map{ [it[0], it[opt_idx]] }.map{ it[1] ? it : [it[0], OPT_FILE_VALUE] }
}

def expand_path( short_path ) {
    return file( short_path ).toAbsolutePath()
}

def get_size_in_gb( files ) {
    if ( files instanceof List ) {
        return (files.sum{ f -> f.size() / 1073741824 }).GB
    }
    def fs = files.size() / 1073741824
    return (fs < 1 ? 1 : fs).GB
}

def prevent_sci_notation ( float_number ) {
    return String.format("%f", float_number)
}

def extract_extension ( f ) {
    return "$f".tokenize(".").subList(1, "$f".tokenize(".").size()).join(".")
}

def copy_and_rename ( fl, prefix, overwrite ) {
    def ext = extract_extension(fl)
    if ( !file("${file(fl).getParent()}/${prefix}.${ext}").exists() || overwrite == "true" )
        file(fl).mklink("${file(fl).getParent()}/${prefix}.${ext}", overwrite: true)
    return file("${file(fl).getParent()}/${prefix}.${ext}")
}

def uniformize_naming ( files_channel, prefix, overwrite ) {
    return files_channel.map{ it ->
        [it[0]] + it.subList(1, it.size()).collect{ i ->
            copy_and_rename(i, "${i.simpleName.split("__")[0]}__$prefix", overwrite)
        }
    }
}

def replace_naming_to_underscore ( files_channel, prefix, overwrite ) {
    return files_channel.map{ it ->
        [it[0]] + it.subList(1, it.size()).collect{ i ->
            def suffix = i.simpleName().tokenize("_")[-1]
            copy_and_rename(i, "${prefix}_${suffix}", overwrite)
        }
    }
}

def sort_by_name ( channel, reg_list ) {
    return channel.map{ [it[0]] + it[1].sort{ f -> f_token = file("$f").simpleName(); reg_list.find{ pt -> f_token ==~ pt } } }
}

def sort_by_extension ( channel, ext_list ) {
    return channel.map{ [it[0]] + it[1].sort{ f -> f_token = "$f".tokenize('.'); ext_list.indexOf(f_token.subList(1, f_token.size()).join('.')) } }
}

def swap_configurations ( base_config, new_config ) {
    if ( new_config && !new_config.empty() )
        return new_config.name
    return base_config
}

def sort_as_with_name ( channel, sorting_channel ) {
    channel.map{ [it[0], it.subList(1, it.size())] }.join(
        sorting_channel.map{ [it[0], it.subList(1, it.size())] }
    ).map{
        [it[0]] + it[1].sort{ f -> f_token = file("$f").getSimpleName(); it[2].find{ pt -> f_token ==~ pt } }
    }
}

def merge_repetitions ( channel, keep_rep_key ) {
    c = channel.map{ it ->
        def sub_rep = it[0].split("_");
        [sub_rep[0], sub_rep[1..<sub_rep.size()].join("_")] + it.subList(1, it.size())
    }.groupTuple()
    c = c.map{
        [it[0]] + it.subList(keep_rep_key ? 1 : 2, it.size()).collect{
            s -> s.withIndex().collect{
                o, i -> [o: o, i: i]
            }.sort{
                a, b -> it[1][a.i] <=> it[1][b.i]
            }.collect{ e -> e.o }
        }
    }

    return c
}

def interleave ( l1, l2 ) {
    def result = [l1, l2].transpose()
    return ( result += (l1 - result*.get(0)) ?: (l2 - result*.get(1)) ).flatten()
}

def merge_channels_non_blocking ( c1, c2 ) {
    c3 = c1.map{ [it[0], it.subList(1, it.size())] }.join(c2.map{ [it[0], it.subList(1, it.size())] })
    return c3.map{ [it[0]] + it.subList(1, 3).transpose() }
}
