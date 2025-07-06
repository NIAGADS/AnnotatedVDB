"""! @brief Annotated VDB VCF Loader"""

#!pylint: disable=invalid-name

##
# @package loaders
# @file vcf_variant_loader.py
#
# @brief  Annotated VDB VCF Loader
#
# @section vcf_variant_loader Description
# provides functions & class for managing variant loads from a VCF entry
#
# @section todo_vcf_variant_loader TODO
# - create __verify_current_variant() to make sure _current_variant is set before trying to access
# - extract common elements to parent class and make VEPLoader a child
#
# @section libraries_vcf_variant_loader Libraries/Modules
# - json - parse/serialize JSON
# - copy - for deep copy functionality
# - [GenomicsDBData.Util.utils](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/utils.py)
#   + provides variety of wrappers for standard file, string, and logging operations
# - [GenomicsDBData.Util.list_utils](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/list_utils.py)
#   + provides variety list and set operations
# - [GenomicsDBData.Util.postgres_dbi](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/postgres_dbi.py)
#   + wrapper for connecting to database & handling PG errors
# - [AnnotatedVDB.Util.parsers](https://github.com/NIAGADS/AnnotatedVDB/blob/master/Util/lib/python/parsers)
#   + parsers for VCF and VEP JSON
# - [AnnotatedVDB.Util.loaders](https://github.com/NIAGADS/AnnotatedVDB/blob/master/Util/lib/python/loaders)
#   + parent class
# - [AnnotatedVDB.Util.variant_annotator](https://github.com/NIAGADS/AnnotatedVDB/blob/master/Util/lib/python/variant_annotator.py)
#   + generate standard variant annotations from position and allele information
# @section author_vcf_variant_loader Author(s)
# - Created by Emily Greenfest-Allen (fossilfriend) 2022

import json
import re
import subprocess
import traceback

from logging import StreamHandler

from copy import deepcopy
from io import StringIO

from GenomicsDBData.Util.utils import (
    xstr,
    warning,
    print_dict,
    to_numeric,
    deep_update,
    execute_cmd,
)
from GenomicsDBData.Util.list_utils import qw, is_subset, is_equivalent_list
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.parsers import VcfEntryParser, CONSEQUENCE_TYPES
from AnnotatedVDB.Util.variant_annotator import VariantAnnotator
from AnnotatedVDB.Util.loaders import (
    VariantLoader,
    REQUIRED_COPY_FIELDS,
    JSONB_UPDATE_FIELDS,
    BOOLEAN_FIELDS,
)
from AnnotatedVDB.Util.database import VARIANT_ID_TYPES

NBSP = " "  # for multi-line sql


class VCFVariantLoader(VariantLoader):
    """! functions for loading variants from a VCF file"""

    def __init__(self, datasource, verbose=False, debug=False):
        """! VCFVariantLoader base class initializer

        @param datasource          datasource description
        @param logFileHandler      log file handler
        @param verbose             flag for verbose output
        @param debug               flag for debug output

        @returns                   An instance of the VCFVariantLoader class with initialized counters and copy buffer
        """

        self.__update_value_generator = self.__default_update_values
        self.__update_existing = False
        self.__update_fields = None
        self.__vcf_header_fields = None
        self.__chromosome = None
        self.__requireSequenceValidation = True
        self.__structural_variant = False
        self.__pk_map_file = None
        self.__adsp_release = None
        super(VCFVariantLoader, self).__init__(datasource, verbose, debug)
        self.logger.info(type(self).__name__ + " initialized")

    def set_adsp_release(self, rls: str):
        self.__adsp_release = rls.lower()

    def set_chromosome(self, chrom):
        self.__chromosome = chrom

    def structural_variants(self):
        self.__structural_variant = True

    def set_pk_map_file(self, fn):
        self.__pk_map_file = fn

    def require_sequence_validation(self, flag: bool):
        self.__requireSequenceValidation = flag

    def chromosome(self):
        return self.__chromosome

    def set_vcf_header_fields(self, fields):
        self.__vcf_header_fields = fields

    def vcf_header_fields(self):
        return self.__vcf_header_fields

    def set_update_existing(self, updateExisting):
        self.__update_existing = updateExisting

    def update_existing(self):
        return self.__update_existing

    def initialize_copy_sql(self, copyFields=None):
        fields = REQUIRED_COPY_FIELDS  # "chromosome", "record_primary_key", "position", "metaseq_id", "bin_index", "row_algorithm_id"]
        fields.extend(
            [
                "ref_snp_id",
                "is_multi_allelic",
                "display_attributes",
                "allele_frequencies",
            ]
        )

        if self.__structural_variant:
            fields.extend(["is_structural_variant", "vcf_entry", "adsp_qc"])

        if copyFields:
            fields.extend(copyFields)

        if self._debug:
            self.logger.debug("Copy fields: %s", fields)

        super(VCFVariantLoader, self).initialize_copy_sql(copyFields=fields)

    def set_update_fields(self, fields):
        self.__update_fields = deepcopy(fields)

    def set_update_value_generator(self, func):
        self.__update_value_generator = func

    def generate_update_values(self, entry, flag=None):
        return self.__update_value_generator(self, entry, flag)

    def __default_update_values(self, entry, flags):
        # flags should be {"metaseq_id":<value>}
        # check if duplicates
        if self.update_existing():
            raise NotImplementedError(
                "Cannot update existing variants with default generator. See update_from_qc_vcf_file.py for example of custom generator"
            )

    def build_update_sql(self):
        """generate update sql"""
        fields = self.__update_fields
        if fields is None:
            raise ValueError(
                "Need to specify update fields using loader.set_update_fields(<array>) to do an update"
            )

        updateString = ""
        for f in fields:
            if f in JSONB_UPDATE_FIELDS:
                # updateString += ' '.join((f,  '=', 'COALESCE(v.' + f, ",'{}'::jsonb)", '||', 'd.' + f + '::jsonb')) + ',' + NBSP
                updateString += (
                    " ".join(
                        (
                            f,
                            "=",
                            "jsonb_merge(COALESCE(v." + f,
                            ",'{}'::jsonb),",
                            "d." + f + "::jsonb)",
                        )
                    )
                    + ","
                    + NBSP
                )  # e.g. v.gwas_flags = jsonb_merge(v.gwas_flags, d.gwas_flags::jsonb)
                # e.g. v.gwas_flags = v.gwas_flags || d.gwas_flags::jsonb
            elif f in BOOLEAN_FIELDS:
                updateString += " ".join((f, "=", "d." + f + "::BOOLEAN")) + "," + NBSP
            else:
                updateString += " ".join((f, "=", "d." + f)) + "," + NBSP

        updateString = updateString.rstrip(", ") + NBSP

        if (
            self.is_adsp() and "is_adsp_variant" not in fields
        ):  # specified as datasource
            fields.append("is_adsp_variant")
            updateString += ", v.is_adsp_variant = d.is_adsp_variant" + NBSP

        fields = ["record_primary_key", "chromosome"] + fields
        schema = "AnnotatedVDB"
        # partition = '_chr' + xstr(self.chromosome()) if self.chromosome() is not None else ''
        sql = (
            "UPDATE "
            + schema
            + ".Variant v SET"
            + NBSP
            + updateString
            + "FROM (VALUES %s) AS d("
            + ",".join(fields)
            + ")"
            + NBSP
            + "WHERE v.chromosome = d.chromosome"
            + NBSP
            + "AND v.record_primary_key = d.record_primary_key"
        )

        self.set_update_sql(sql)
        if self._debug:
            self.logger.debug("Update SQL: %s", self._update_sql)

    def __buffer_update_values(self, entry, flags):
        """save udpate values to value list"""
        # TODO: check for duplicate has to happen in generate_update_values if not accounted for by flags
        # see default
        recordPK, uFlags, uValues = self.generate_update_values(entry, flags)
        isAdspVariant = None

        if uFlags is not None:
            if "update" in uFlags and uFlags["update"] is False:
                if self._debug:
                    self.logger.debug(
                        recordPK + " already has values for update fields; skipping"
                    )
                self.increment_counter("skipped")
                return "SKIPPED"

            if "is_adsp_variant" in "uFlags":
                isAdspVariant = uFlags["is_adsp_variant"]

        if recordPK is None:
            if self._debug:
                self.logger.debug(
                    "Variant %s not in DB, inserting", self._current_variant.id
                )
            return "INSERT"

        values = [recordPK]
        values.append("chr" + self._current_variant.chromosome)

        if self._debug:
            self.logger.debug(
                "Updating variant %s (is ADSP = %s) with values: %s",
                recordPK,
                xstr(isAdspVariant),
                uValues,
            )

        for f in self.__update_fields:
            cValue = uValues[f]

            if f in JSONB_UPDATE_FIELDS:
                values.append(xstr(cValue))
            else:
                values.append(uValues[f] if uValues[f] != "NULL" else None)

        if self.is_adsp() and "is_adsp_variant" not in uValues:
            values.append(True)

        if self._debug:
            self.logger.debug("Buffered Values: %s", values)

        self._update_buffer.append(tuple(values))
        self.increment_counter("update")
        return "UPDATE"

    def __add_adsp_update_statement(self, recordPK, chromosome):
        """add update is_adsp_variant to update buffer"""
        chrm = "chr" + xstr(chromosome)
        values = (recordPK, chrm)
        self._update_buffer.append(values)
        # updateStr = "UPDATE AnnotatedVDB.Variant SET " \
        #    + "is_adsp_variant = true WHERE record_primary_key = '"  + recordPK + "' " \
        #    + "AND chromosome = '" + chrm + "'" # incl chromosome so it uses an index
        # self._update_buffer.write(updateStr + ';')

    def __generate_primary_key(
        self, chrm, pos, ref, alt, externalId, requireValidation=True
    ):
        """
        wrapper for generate primary key to catch invalid indels
        """

        try:
            annotator = VariantAnnotator(ref, alt, chrm, pos)
            # rv = True if self.__requireSequenceValidation == True else requireValidation
            recordPK = self._pk_generator.generate_primary_key(
                annotator.get_metaseq_id(),
                externalId,
                requireValidation=requireValidation,
            )
        except ValueError as err:
            # sequence mismatch for long indel; check to see if possibly an insertion
            try:
                annotator = VariantAnnotator(alt, ref, chrm, pos)
                recordPK = self._pk_generator.generate_primary_key(
                    annotator.get_metaseq_id(), externalId, requireValidation=True
                )
                self.logger.warning(f"switching alleles: " + str(err))

            except:
                if self.debug:
                    self.logger.debug("Final lookup")
                annotator = VariantAnnotator(ref, alt, chrm, pos)
                recordPK = self._pk_generator.generate_primary_key(
                    annotator.get_metaseq_id(), externalId, requireValidation=False
                )
        return annotator, recordPK

    def __generate_sv_primary_key(self, chrom, start, end, svType):
        recordPK = self._pk_generator.generate_sv_primary_key(chrom, start, end, svType)
        return recordPK

    def __is_duplicate_from_file(self, pk):
        """
        Search for a pk in self.__pk_map_file.
        Returns the matching line or None if not found.
        """
        if not self.__pk_map_file:
            raise ValueError("Primary key map file not set.")

        cmd = ["grep", "-m", "1", pk, self.__pk_map_file]
        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            if result.returncode == 0 and result.stdout:
                self.logger.debug(f"{result.stdout}")
                return True
            else:
                return False
        except Exception as e:
            raise e

    def __parse_structural_variant(self, vcfEntry: VcfEntryParser, flags):
        if self._debug and self._verbose:
            self.logger.debug("Current Variant: %s", print_dict(self._current_variant))
            self.logger.debug("Update/Load flags: %s", flags)

        chrom = vcfEntry.get("chrom")
        start = vcfEntry.get("pos")
        ref = vcfEntry.get("ref")
        end = vcfEntry.get_info("END")
        svType = vcfEntry.get_info("SVTYPE")

        if ref != "N":
            self.increment_counter("skipped")
            if self._debug:
                self.logger.debug(f"Skipping non-SV {vcfEntry.get("id")}")
            return None

        recordPK = self.__generate_sv_primary_key(chrom, start, end, svType)

        if self._skip_existing:
            if self.__is_duplicate_from_file(recordPK):
                if self._debug:
                    self.logger.debug(f"Skipping duplicate {recordPK}")
                self.increment_counter("skipped")
                return None

            if self.is_duplicate(recordPK, returnMatch=True):
                if self._debug:
                    self.logger.debug(f"Skipping duplicate {recordPK}")
                self.increment_counter("skipped")
                return None

        binIndex = self._bin_indexer.find_bin_index(chrom, start, end)

        qcScores = deepcopy(vcfEntry.get("info"))
        adspQC = {
            self.__adsp_release: {
                "info": qcScores,
                "filter": vcfEntry.get("filter"),
                "qual": vcfEntry.get("qual"),
                "format": vcfEntry.get("format"),
            }
        }

        svSize = vcfEntry.get_info("SVSIZE")
        if svSize is None:
            # alt = "<INS:SVSIZE=121+:AGGREGATED>"
            match = re.search(r"SVSIZE=([0-9+]+)", self._current_variant.alt_alleles[0])
            if match:
                svSize = match.group(1).replace("+", "")

        svModel = vcfEntry.get_info("SVMODEL")
        if svModel is None:
            # alt = "<INS:SVSIZE=121+:AGGREGATED>"
            size, model = self._current_variant.alt_alleles[0].split(":")
            svSize = model.replace(">", "")

        display = {
            "location_start": int(start),
            "location_end": int(end),
            "variant_class_abbrev": svType,
            "sv_size": xstr(svSize, nullStr="NULL"),
            "sv_model": svModel,
        }

        copyValues = [
            "chr" + xstr(self._current_variant.chromosome),
            recordPK,
            xstr(self._current_variant.position),
            recordPK,
            binIndex,
            xstr(self._alg_invocation_id),
            "NULL",  # refsnp
            "NULL",  # is_multiallelic
            xstr(display, nullStr="NULL"),
            "NULL",  # allele frequency
            xstr(True),  # is_structural_variant
            # passing NULL for these for now b/c QC causes duplicates
            # xstr(vcfEntry.get_entry(), nullStr="NULL"),  # vcf_entry
            # xstr(adspQC, nullStr="NULL"),  # adsp_qc
            "NULL",
            "NULL",
        ]

        if self.is_adsp():
            copyValues.append(xstr(True))  # is_adsp_variant

        if self._debug:
            self.logger.debug("COPY values: %s", copyValues)

        self.add_copy_str("#".join(copyValues))
        self.increment_counter("variant")

        primaryKeyMapping = [{"primary_key": recordPK, "bin_index": binIndex}]
        return {self._current_variant.id: primaryKeyMapping}
        # TODO raise error on duplicate PKs

    def __parse_alt_alleles(self, vcfEntry, flags):
        """! iterate over alleles and calculate values to load
        add to copy buffer
        @param vcfEntry             the VCF entry for the current variant
        @param flags                info passed for script that may qualify update or load
        """

        if self._debug and self._verbose:
            self.logger.debug("Current Variant: %s", print_dict(self._current_variant))
            self.logger.debug("Update/Load flags: %s", flags)

        variant = self._current_variant  # just for ease/simplicity
        variantExternalId = (
            variant.ref_snp_id if hasattr(variant, "ref_snp_id") else None
        )

        primaryKeyMapping = []

        for alt in variant.alt_alleles:
            if alt == ".":
                self.logger.warning(
                    "Skipping variant " + variant.id + "; no alt allele (alt = .)"
                )
                self.increment_counter("skipped")
                continue

            annotator, recordPK = self.__generate_primary_key(
                variant.chromosome,
                variant.position,
                variant.ref_allele,
                alt,
                variantExternalId,
            )

            matchedVariant = None
            if self.skip_existing():
                matchedVariant = self.is_duplicate(
                    annotator.get_metaseq_id(), returnMatch=True
                )
                if matchedVariant:
                    primaryKeyMapping += matchedVariant
                    if self._log_skips:
                        self.logger.info(
                            "Duplicate found %s: %s",
                            annotator.get_metaseq_id(),
                            matchedVariant,
                        )
                    self.increment_counter("skipped")

            if matchedVariant is None:
                status = None
                if flags is None:
                    flags = {"metaseq_id": annotator.get_metaseq_id()}
                if self.update_existing():
                    status = self.__buffer_update_values(vcfEntry, flags)
                    if status != "INSERT":
                        continue  # skipped or updated

                # FIXME: don't update ADSP status here/ remove this from the code
                if self.is_adsp():
                    if self.is_duplicate(recordPK):
                        self.__add_adsp_update_statement(
                            recordPK, self._current_variant.chromosome
                        )
                        self.increment_counter("update")
                        continue

                normRef, normAlt = annotator.get_normalized_alleles()
                binIndex = self._bin_indexer.find_bin_index(
                    variant.chromosome,
                    variant.position,
                    vcfEntry.infer_variant_end_location(alt),
                )

                alleleFreq = vcfEntry.get_frequencies(alt)

                additionalValues = None
                if self.__update_fields is not None:
                    pkPlaceholder, uFlags, additionalValues = (
                        self.generate_update_values(vcfEntry, flags)
                    )

                # TODO iterate over copy fields and add as each one is hit so order is preserved
                copyValues = [
                    "chr" + xstr(variant.chromosome),
                    recordPK,
                    xstr(variant.position),
                    annotator.get_metaseq_id(),
                    binIndex,
                    xstr(self._alg_invocation_id),
                    xstr(variant.ref_snp_id, nullStr="NULL"),
                    xstr(variant.is_multi_allelic, falseAsNull=True, nullStr="NULL"),
                    xstr(
                        annotator.get_display_attributes(variant.rs_position),
                        nullStr="NULL",
                    ),
                    xstr(alleleFreq, nullStr="NULL"),
                ]

                if additionalValues:
                    for f in self.__update_fields:
                        copyValues.append(xstr(additionalValues[f], nullStr="NULL"))

                if self.is_adsp():
                    copyValues.append(xstr(True))  # is_adsp_variant

                if self._debug:
                    self.logger.debug(
                        "Display Attributes: "
                        + print_dict(
                            annotator.get_display_attributes(variant.rs_position)
                        )
                    )
                    self.logger.debug("COPY values: %s", copyValues)

                self.add_copy_str("#".join(copyValues))
                self.increment_counter("variant")

                primaryKeyMapping += [{"primary_key": recordPK, "bin_index": binIndex}]

        return {self._current_variant.id: primaryKeyMapping}

    def parse_variant(self, line, flags=None):
        """! parse & load single line from file
        @param line             the VCF entry in string form
        @param flags            flags that may be used to modify update or copy strings
        @returns copy string for db load
        """

        if self.resume_load() is False and self._resume_after_variant is None:
            raise ValueError(
                "Must set VariantLoader resume_afer_variant if resuming load"
            )

        variantPrimaryKeyMapping = None

        try:
            self.increment_counter("line")
            entry = (
                VcfEntryParser(line, self.vcf_header_fields())
                if isinstance(line, str)
                else line
            )  # assume vcfEntry was passed

            if not self.resume_load():
                self._update_resume_status(entry.get("id"))
                return None

            # otherwise proceed with the parsing
            entry.update_chromosome(self._chromosome_map)

            # extract identifying variant info for frequent reference & logging in the loader calling script
            self._current_variant = entry.get_variant(
                dbSNP=self.is_dbsnp(), namespace=True
            )

            if self.__structural_variant:
                if entry.get("filter") == "PASS":  # for now
                    variantPrimaryKeyMapping = self.__parse_structural_variant(
                        entry, flags
                    )
                else:
                    if self._debug:
                        self.logger.debug(f"Skipping non-PASS record {entry.get("id")}")
                    self.increment_counter("skipped")
                    return None

                if variantPrimaryKeyMapping is None:  # duplicate
                    return None

            else:
                # iterate over alleles
                variantPrimaryKeyMapping = self.__parse_alt_alleles(entry, flags)

            if self._debug:
                self.logger.debug("PKM: %s", variantPrimaryKeyMapping)

            if variantPrimaryKeyMapping is None:
                raise ValueError(entry)

            return variantPrimaryKeyMapping

        except Exception as err:
            raise err
