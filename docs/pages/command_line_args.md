# Command line arguments

:::{container} llm-review-note
**LLM-authored draft — review required.** This page was drafted from the 2025 installation/use manual and checked against the available command-line options. Review the operational guidance before release.
:::

The installation manual describes IdrAgra as a standalone command-line program tested on Windows and Linux, with no additional runtime libraries required by the distributed build from version 2 onward. Building this repository from source still requires the compiler and tools described by the repository build setup.

Run IdrAgra from a terminal whose working directory is the project folder. The executable filename depends on the distribution: the manual uses `idragra.exe`, while this repository's Makefile also creates `release\idragra_latest.exe`. With no arguments, it reads `idragra_parameters.txt` from the working directory:

```powershell
C:\path_to_idragra\idragra.exe
```

To select another parameter file, use `-f` (or `-filename`):

```powershell
C:\path_to_idragra\idragra.exe -f .\scenarios\dry_year.txt
```

Relative paths inside the selected parameter file are still resolved from the process working directory, not from the parameter file's own folder.

## Supported options

| Short form | Long form | Current behavior |
|---|---|---|
| `-h` | `-help` | Print the option list and stop. |
| `-d` | `-default` | Print all compiled default parameter values and stop. |
| `-p` | `-preview` | Print the effective settings after the parameter file has been read, then continue the run. |
| `-v` | `-verbose` | Print additional diagnostic and progress messages. It does not enable result files. |
| `-s` | `-summary` | Write only irrigation maps from the normal periodic and annual map groups. Independently enabled yield, debug, and cell outputs are still written. |
| `-f <file>` | `-filename <file>` | Read settings from `<file>` instead of `idragra_parameters.txt`. |

Options are case-insensitive. An unrecognized option that begins with `-` prints the help text and stops.

:::{container} manual-code-divergence
**Manual/code divergence to review.** The installation manual shows a parameter filename as a bare argument, but `idragra.exe my_parameters.txt` still reads the default file. Use `-f my_parameters.txt`.
:::

:::{container} manual-code-divergence
**Manual/code divergence to review.** The manuals list `-t`, `-teta`, and `-theta`, but these options are not supported. Use {ref}`FinalThetaFlag <parameter-finaltheta>`, {ref}`FinalConditionPath <parameter-finalconditionpath>`, and {ref}`FinalCondition <parameter-finalcondition>` to write final soil-moisture grids.
:::

## Output-folder prompt

When {ref}`OutputPath <parameter-outputpath>` points to an existing directory, IdrAgra reports that it will update the directory and waits for Enter. Existing files with names reused by the run are overwritten; unrelated files are left in place. Automated runs should therefore prepare a new output folder or provide input to the process.
