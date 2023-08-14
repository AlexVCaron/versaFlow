#!/usr/bin/env nextflow

package groovyx.gpars.dataflow

nextflow.enable.dsl=2

def group_channel_rep ( chan ) {
    return chan.groupTuple().map{
        [it[0]] + it[1..-1].inject((1..it.size()).collect{ [] }) { sub, rep ->
            sub.eachWithIndex { its, i -> its.add(rep[i]) } ; return sub
        }
    }
}

def group_subject_reps ( dwi_channel, metadata_channel ) {
    return [group_channel_rep(dwi_channel), group_channel_rep(metadata_channel)]
}

def replace_dwi_file ( base_channel, dwi_channel ) {
    return dwi_channel.join(base_channel.map{ [it[0]] + it[2..-1] })
}

OPT_FILE_VALUE = ""
OPT_CHANNEL = null

def opt_channel () {
    return OPT_CHANNEL
}

def is_data ( opt_channel ) {
    return ( opt_channel && !( opt_channel instanceof DataflowVariable ) ) || ( opt_channel instanceof DataflowVariable && opt_channel.value )
}

def join_optional ( base_channel, opt_channel ) {
    if ( is_data(opt_channel) )
        return base_channel.join(opt_channel)
    else
        return base_channel.map { it + [OPT_FILE_VALUE] }
}

def map_optional ( base_channel, opt_idx ) {
    return base_channel.map{ [it[0], it[opt_idx]] }.map{ it[1] ? it : [it[0], OPT_FILE_VALUE] }
}

def expand_path( short_path ) {
    return file( short_path ).toAbsolutePath()
}

def prevent_sci_notation ( float_number ) {
    return String.format("%f", float_number)
}

def extract_extension ( f ) {
    return "$f".tokenize(".")[1..-1].join(".")
}

def copy_and_rename ( fl, prefix, overwrite, copy ) {
    def ext = extract_extension(fl)
    if ( !file("${prefix}.${ext}").exists() || overwrite == "true" )
        if ( copy == "true" )
            file(fl).copyTo("${prefix}.${ext}")
        else
            file(fl).mklink("${prefix}.${ext}", overwrite: true)
    return file("${prefix}.${ext}")
}

def uniformize_naming ( files_channel, prefix, overwrite, copy ) {
    return files_channel.map{ it ->
        [it[0]] + it[1..-1].collect{ i ->
            i == "" ? i : copy_and_rename(i, "${i.simpleName.split("__")[0]}__$prefix", overwrite, copy)
        }
    }
}

def rename_according_to ( file_channel, ref_channel, suffix, overwrite ) {
    return file_channel.join(ref_channel).map{ it ->
        [it[0]] + it[1..-1].collect{ i ->
            copy_and_rename(i, "${it[-1].simpleName.split("__")[0]}__$suffix", overwrite)
        }
    }
}

def replace_naming_to_underscore ( files_channel, prefix, overwrite ) {
    return files_channel.map{ it ->
        [it[0]] + it[1..-1].collect{ i ->
            def suffix = i.simpleName().tokenize("_")[-1]
            copy_and_rename(i, "${prefix}_${suffix}", overwrite)
        }
    }
}

def sort_by_name ( channel, reg_list ) {
    return channel.map{ [it[0]] + it[1].sort{ f -> f_token = file("$f").simpleName(); reg_list.find{ pt -> f_token ==~ pt } } }
}

def sort_by_extension ( channel, ext_list ) {
    return channel.map{ [it[0]] + it[1].sort{ f -> f_token = "$f".tokenize('.'); ext_list.indexOf(f_token[1..-1].join('.')) } }
}

def swap_configurations ( base_config, new_config ) {
    if ( new_config && !new_config.empty() )
        return new_config.name
    return base_config
}

def sort_as_with_name ( channel, sorting_channel ) {
    channel.map{ [it[0], it[1..-1]] }.join(
        sorting_channel.map{ [it[0], it[1..-1]] }
    ).map{
        [it[0]] + it[1].sort{ f -> f_token = file("$f").getSimpleName(); it[2].find{ pt -> f_token ==~ pt } }
    }
}

def merge_repetitions ( channel, keep_rep_key ) {
    def c = channel.map{ it ->
        def sub_rep = it[0].split("_");
        [sub_rep[0], sub_rep[1..<sub_rep.size()].join("_")] + it[1..-1]
    }.groupTuple()
    c = c.map{
        [it[0]] + it[(keep_rep_key ? 1 : 2)..-1].collect{
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
    def c3 = c1.map{ [it[0], it[1..-1]] }.join(c2.map{ [it[0], it[1..-1]] })
    return c3.map{ [it[0]] + it[1..2].transpose() }
}

def remove_alg_suffixes ( f ) {
    def name = f.split("__")
    if (name.size() == 1) return f
    else return [name[0], extract_extension(f)].join(".")
}

def add_suffix ( f, suffix ) {
    return [f.tokenize(".")[0] + suffix, extract_extension(f)].join(".")
}

def rename( channel, name ) {
    return channel.map{
        [it[0], copy_and_rename(it[1], "${it[0]}_$name", false, true)]
    }
}

def fill_missing_datapoints( data_channel, id_channel, filter_index, fill_tuple ) {
    return id_channel.join(data_channel, remainder: true).map{
        it instanceof ArrayList ? it : [it, null]
    }.map{
        (it[filter_index] == null) ? it[0..filter_index - 1] + fill_tuple : it
    }
}

def filter_datapoints( data_channel, filter_function ) {
    return data_channel.filter{ filter_function(it) }
}

def exclude_missing_datapoints( data_channel, filter_index, missing_value ) {
    return filter_datapoints(data_channel, { it[filter_index] != missing_value })
}

def separate_b0_from_dwi( data_channel ) {
    return [exclude_missing_datapoints(data_channel, 2, ""), filter_datapoints(data_channel, { it[2] == "" })]
}

LIBRARY_ROOT_NAME = "versaFlow"

def get_data_path () {
    def current_dir = file("$projectDir")
    while ( current_dir.name != LIBRARY_ROOT_NAME ) {
        current_dir = current_dir.parent
    }
    return "$current_dir/.data"
}

def get_config_path () {
    def current_dir = file("$projectDir")
    while ( current_dir.name != LIBRARY_ROOT_NAME ) {
        current_dir = current_dir.parent
    }
    return "$current_dir/.config"
}

def is_path_list ( pth ) {
    return ( pth instanceof Path ? pth.getNameCount() : pth.size() ) > 0
}

def collect_paths ( data_channel ) {
    return data_channel.map{
        [it[0], (it.size() > 2 ? it[1..-1] : [it[1]] ).findAll{ it }]
    }
}