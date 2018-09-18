# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""reserves/creates some cmdline options that we need placeholders for"""

from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.common.secmet import Record
from antismash.common import subprocessing

from .results import T2PKS_Results
from .t2pks_analysis import t2pks_analysis
from .html_output import will_handle, generate_sidepanel

NAME = "t2_pks"
SHORT_DESCRIPTION = "type II PKS analysis"

# because this module has a lot of stubbed functions, quiet the pylint warnings
# about unused arguments
# pylint: disable=unused-argument


def get_arguments() -> ModuleArgs:
    """ Constucts T2 PKS module arguments
    """
    args = ModuleArgs("T2 PKS", "t2pks", enabled_by_default=True)
    return args


def check_options(options: ConfigType) -> List[str]:
    """ Check the options of this module for any conflicting or invalid values.

        Arguments:
            options: the options parsed by the main entry point as an
                     antismash.Config object

        Returns:
            a list of strings describing any errors, if they exist
    """
    # because the placeholders are special, raise an error if they're used
    if options.without_fimo:
        raise ValueError("Dummy options can't be enabled")
    return []


def check_prereqs() -> List[str]:
    """ Check the prerequisites.
            hmmscan: domain detection
            blastp: CLF and starter unit analysis
            HMMs: t2pks.hmm

        Returns:
            a list of strings describing any errors, if they occurred
    """
    failure_messages = []
    for binary_name in ['hmmscan', 'blastp']:
        if path.locate_executable(binary_name) is None:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    for hmm in ['t2pks.hmm']:
        hmm = path.get_full_path(__file__, 'data', hmm)
        if path.locate_file(hmm) is None:
            failure_messages.append("Failed to locate file %r" % hmm)
            continue
        for ext in ['.h3f', '.h3i', '.h3m', '.h3p']:
            binary = "%s%s" % (hmm, ext)
            if path.locate_file(binary) is None:
                # regenerate them
                result = subprocessing.run_hmmpress(hmm)
                if not result.successful():
                    failure_messages.append("Failed to hmmpress %s: %s" % (hmm, result.stderr.rstrip()))
                break
    
    for blastdb in ['KSIII', 'AT', 'LIG']:
        for ext in ['.fasta','.phr', '.pin', '.psq']:
            dbfile = path.get_full_path(__file__, 'data', blastdb + ext)
            if path.locate_file(dbfile) is None:
                failure_messages.append("Failed to locate file %r" % dbfile)
                continue

    return failure_messages


def is_enabled(options: ConfigType) -> bool:
    """ No options to check at the moment """
    return not options.minimal or options.t2_pks_enabled


def regenerate_previous_results(previous: Dict[str, Any], record: Record, options: ConfigType) -> T2PKS_Results:
    """ Rebuild the previous run results from a JSON object into this module's
        python results class. If the current options are incompatible with the
        previous results, None should be returned.

        The module result class should inherit from
            antismash.common.module_results.ModuleResults

        This module doesn't ever run, so it doesn't have any results to
        regenerate.

        Arguments:
            previous: the previous results as a dictionary
            record: the Record that was used to generate the previous results
            options: an antismash.Config object
    """
    return T2PKS_Results.from_json(previous, record)


def run_on_record(record: Record, results: Optional[T2PKS_Results], options: ConfigType) -> T2PKS_Results:
    """ Run this module's analysis section on the given record or use the
        previous results.

        Arguments:
            record: the Record instance to analyse
            results: the previous results as generated by regenerate_previous_results()
            options: an antismash.Config object

        Returns:
            this module's results as a subclass of
                antismash.common.module_results.ModuleResults
    """
    if not results:
        results = T2PKS_Results(record.id)
        t2pks_analysis(record, results, options)

    return results
    
    #raise NotImplementedError("Dummy module should never be run")