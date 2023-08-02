import os
import argparse
import pathlib
from typing import Sequence, Union
from input_prep.query_utils import parse_fasta_files, parse_csv_files, clean_and_validate_queries, detect_duplicate_queries

class FileArgumentParser(argparse.ArgumentParser):
    """Overwrites default ArgumentParser to better handle flag files."""

    def convert_arg_line_to_args(self, arg_line: str) -> Sequence[str]:
        """ Read from files where each line contains a flag and its value, e.g.
        '--flag value'. Also safely ignores comments denotes with '#' and
        empty lines.
        """

        # Remove any comments from the line
        arg_line = arg_line.split('#')[0]

        # Escape if the line is empty
        if not arg_line:
            return None

        # Separate flag and values
        split_line = arg_line.strip().split(' ')

        # If there is actually a value, return the flag value pair
        if len(split_line) > 1:
            return [split_line[0], ' '.join(split_line[1:])]
        # Return just flag if there is no value
        else:
            return split_line

class QueryManager(object):
    """Manager that will parse, validate, and store queries. """

    def __init__(self,
            input_dir: str = None,
            sequences: Sequence[Union[str, Sequence[str]]] = [],
            min_length: int = 16,
            max_length: int = 2500,
            max_multimer_length: int = 2500) -> None:

        self.sequences = sequences
        
        self.min_length = min_length
        self.max_length = max_length
        self.max_multimer_length = max_multimer_length

        self.files = {}
        self.others = []
        self.queries = []

        # Detect .fasta and .csv files from the input directory.
        if input_dir not in ['', None]:
            onlyfiles = [f for f in os.listdir(input_dir) if os.path.isfile(
                         os.path.join(input_dir, f))]

            for filename in onlyfiles:
                extension = filename.split('.')[-1]
                if extension in ['fasta', 'csv']:
                    if extension not in self.files:
                        self.files[extension] = []
                    self.files[extension].append(os.path.join(input_dir, filename))
                else:
                    self.others.append(os.path.join(input_dir, filename))
                
        if len(self.files) == 0 and self.sequences == []:
            raise ValueError(
                f'No input files (.fasta or .csv) detected in '
                '{input_dir} and no sequences provided.')

        if len(self.files) != 0:
            filenames = [pathlib.Path(p).stem for l in self.files.values() for p in l]
            if len(filenames) != len(set(filenames)):
                raise ValueError('All input files must have a unique basename.')
        
    def parse_files(self) -> None:

        # For each detected filetype parse queries
        for extension in self.files:
            if extension == 'fasta':
                queries = parse_fasta_files(
                    files=self.files['fasta'])
            else:
                queries = parse_csv_files(
                    files=self.files['csv'])

            # Validate queries by checking sequence composition and lengths
            queries = clean_and_validate_queries(
                input_queries=queries,
                min_length=self.min_length,
                max_length=self.max_length,
                max_multimer_length=self.max_multimer_length)

            # Add queries to overall query list.
            self.queries += queries
                            
        # Remove duplicate queries to reduce redundancy
        self.queries = detect_duplicate_queries(
            query_list=self.queries)

    def parse_sequences(self) -> None:

        queries = []
        
        for sequence in self.sequences:
            queries.append( ('_INPUT_', sequence) )

        #print(queries)
        if queries != []:
            queries = clean_and_validate_queries(
                input_queries=queries,
                min_length=self.min_length,
                max_length=self.max_length,
                max_multimer_length=self.max_multimer_length)

        # Add queries to overall query list.
        self.queries += queries

        # Remove duplicate queries to reduce redundancy.
        self.queries = detect_duplicate_queries(
            query_list=self.queries)