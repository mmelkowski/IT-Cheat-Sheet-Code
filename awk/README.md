# AWK

The AWK language is a data-driven scripting language consisting of a set of actions to be taken against streams of textual data.

TLDR: It's great but it's complicated and it do things line by line

## Examples of quick command

### Select row on conditions

#### General

```bash
awk '{ if(condition) { print }}' inputfile > outputfile
```

#### Multiple conditions

```bash
awk '{ if ((condition) && (condition) && (condition)) { print } }' inputfile > outputfile
```

#### Example

```bash
awk '{ if($4 >= 5000) { print }}' inputfile > outputfile
```

### Extract column on sep

#### General

```bash
awk -F'separator' '{print $num_column}' inputfile > outputfile
```

#### Example

```bash
awk -F'\t' '{print $2}' inputfile > outputfile
```
