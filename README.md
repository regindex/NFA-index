# p-sort-index

This tool build an index for a sorted NFA that, support count, locate and membership query. The index implemented is a refined version of [1]

### Requirements

* A modern C++17 compiler such as `g++` version 8.4 or higher.
* The cmake build system, version \geq 3.10.
* A Linux or MacOS 64-bit operating system.
* The sdsl-lite library
* The c++ Boost library, version $\geq$ 3.10

### Download

To clone the repository, run:

```console
git clone http://github.com/regindex/p-sort-index
```

### Compile

You can compile using cmake and make:

```console
mkdir build
cd build
cmake ..
make
```

### Input

This tool accepts a co-lexicographically ordered automaton in .graphml format. The source state must be the smallest state in the first chain.  States have three fields: a boolean indicating whether the state is final, the state's chain ID, and the state's position within that chain. Edges have a field storing their label.  A complete example of an ordered NFA is provided in the ./example folder.
```
<node id="1">
    <data key="final">0</data>
    <data key="chains">1</data>
    <data key="positions">1</data>
</node>
...
<edge source="1" target="7">
    <data key="label">b</data>
</edge>
```
The file supports multiple chain partitions within the same graph.  The chains and positions fields specify these partitions, with comma separation used to delineate the different partitions for a given state.
```
<node id="1">
    <data key="final">0</data>
    <data key="chains">1,1</data>
    <data key="positions">1,2</data>
</node>
```

### Usage

```
usage: ./p-sort-index [-h] [-b BUILD_TYPE] [-o OUT_PATH] [-i INPUT_PATH] [-p OPERATION] [-c CHAIN] [-d] [-s SIZE_PATH]

Tool to sort and prune finite automata using the partition refinement algorithm.

Options:
  -h [ --help ]                Help screen
  -b [ --build ] BUILD_TYPE    Build construction type: 0 read a NFA a build 
                               index, 1 load the index from file
  -i [ --input ] INPUT_PATH    Input path of the NFA or Index
  -o [ --output ] OUTPUT_PATH  Output path for the 0 construction type 
                               (default: index.sdsl)
  -p [ --operation ] OPERATION Operation to be performed: count, locate or 
                               membership
  -c [ --chain ] CHAIN         Use Chain partition (CHAIN) to build the Index 
                               (if loaded from NFA, default 1)
  -d [ --debug ]               Debug options to add some debug output
  -s [ --size ] SIZE_PATH      Print size of the index structures in Bytes 
                               inside the SIZE_PATH file
```

When an operation is selected, input is read from the standard input stream, and the result is written to the standard output stream.

### Example

```
./p-sort-index -b 0 -i ../example/sorted-NFA.graphml -o NFA_index.sdsl
./p-sort-index -b 1 -i NFA_index.sdsl -p count
>ab
2
>cd
1
```

## References

- [1] Nicola Cotumaccio, Giovanna D’Agostino, Alberto Policriti, Nicola Prezza, Co-lexicographically ordering automata and regular languages-part i: https://dl.acm.org/doi/full/10.1145/3607471

### Funding

This project has received funding from the European Research Council (ERC) under the European Union’s Horizon Europe research and innovation programme, project REGINDEX, grant agreement No 101039208
