#!/usr/bin/rdmd

import core.stdc.stdlib : exit;

import std.stdio;
import std.getopt;
import std.format;
import std.file;
import std.process;
import std.array;

string compiler = "dmd";
string[] examples = [];

string[] library_files = [
    "linalg/vector.d",
    "linalg/matrix.d",
    "quantum/qc.d",
    "algos/qft.d"
];

void main(string[] args) {
    auto help_info = getopt(args,
        "c|compiler", "specify the compiler to use (default: dmd)", &compiler,
        "e|example", "specify an example or examples to build", &examples
    );

    if (help_info.helpWanted) {
        defaultGetoptPrinter("QSLD build script", help_info.options);
    }

    const string[] LIBRARY_BUILD_COMMAND = [
        compiler, "-lib", "-O", "-of=libqsld.a"
    ] ~ library_files;

    string f = "libqsld.a";
    if (!f.exists) {
        auto compilation_pid = spawnProcess(LIBRARY_BUILD_COMMAND);
        if (wait(compilation_pid) != 0) {
            writeln("Compilation of libqsld.a failed");
            exit(1);
        }
    }

    foreach (example; examples) {
        string bin_name_path = format("./examples/%s", example);
        string src_name_path = format("./examples/%s.d", example);

        string[] compile_example_cmd = [
            compiler, "-O", "-L-L.", "-L=-lqsld",
            "-of=" ~ bin_name_path,
            src_name_path,
        ];

        if (!src_name_path.exists) {
            writeln("The example source file with the name ", src_name_path, " does not exist...skipping");
            continue;
        }

        if (!bin_name_path.exists) {
            auto compilation_pid = spawnProcess(compile_example_cmd);
            if (wait(compilation_pid) != 0) {
                writeln("Compilation of example with name ", bin_name_path, " failed");
                exit(1);
            }
        } else {
            writeln("The example binary with name ", bin_name_path, " already exists...skipping");
        }
    }
}
