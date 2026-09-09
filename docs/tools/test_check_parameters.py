"""Regression tests for the parameter documentation checker."""

from __future__ import annotations

import contextlib
import io
import unittest

import check_parameters as checker


class ParameterCheckerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        parsed, cls.parser_duplicates = checker.parse_parser(checker.DEFAULT_PARSER)
        declarations = checker.parse_declarations(checker.DEFAULT_DECLARATIONS)
        cls.parameters = checker.enrich_parameters(parsed, declarations)
        cls.by_name = {parameter.name: parameter for parameter in cls.parameters}

    def test_outer_parameter_cases_are_extracted_without_nested_enum_cases(self) -> None:
        self.assertFalse(self.parser_duplicates)
        self.assertIn("monthlyflag", self.by_name)
        self.assertNotIn("monthly", self.by_name)
        self.assertIn("weekday", self.by_name)
        self.assertNotIn("monday", self.by_name)

    def test_regular_scalar_metadata_is_resolved(self) -> None:
        parameter = self.by_name["meteostatweightnum"]
        self.assertEqual(parameter.type_name, "Integer")
        self.assertEqual(parameter.default, "1")

    def test_array_element_defaults_follow_fortran_column_major_order(self) -> None:
        self.assertEqual(self.by_name["startdate"].default, "10")
        self.assertEqual(self.by_name["09q_eva"].default, "8.026400")
        self.assertEqual(self.by_name["01q_trasp"].default, "0.472116")

    def test_uninitialized_declaration_has_no_default(self) -> None:
        parameter = self.by_name["finalthetaflag"]
        self.assertEqual(parameter.type_name, "Boolean")
        self.assertEqual(parameter.default, "none")

    def test_documented_metadata_is_read_from_compact_entries(self) -> None:
        documented, duplicates = checker.parse_documentation(checker.DEFAULT_DOCUMENTATION)
        by_name = {parameter.name: parameter for parameter in documented}
        self.assertFalse(duplicates)
        self.assertEqual(by_name["outputpath"].type_name, "String")
        self.assertEqual(by_name["meteostatweightnum"].default, "`1`")
        self.assertEqual(by_name["lim_prec"].type_name, "Real")
        self.assertEqual(by_name["finalthetaflag"].default, "none")

    def test_translated_parameters_use_their_user_facing_formats(self) -> None:
        self.assertEqual(
            (self.by_name["monthlyflag"].type_name, self.by_name["monthlyflag"].default),
            ("String", "monthly"),
        )
        self.assertEqual(
            (self.by_name["weekday"].type_name, self.by_name["weekday"].default),
            ("String or Integer", "monday"),
        )
        self.assertEqual(
            (self.by_name["simulatedsoiluses"].type_name, self.by_name["simulatedsoiluses"].default),
            ("Integer array", "none"),
        )

    def test_reference_contains_every_parameter_assigned_by_the_demo(self) -> None:
        demo_path = checker.REPOSITORY_ROOT / "demo" / "idragra_parameters.txt"
        demo_names = {
            line.split("=", 1)[0].strip().casefold()
            for line in checker.read_text(demo_path).splitlines()
            if "=" in line and not line.lstrip().startswith("#")
        }
        documented, _ = checker.parse_documentation(checker.DEFAULT_DOCUMENTATION)
        self.assertEqual({item.name for item in documented}, demo_names)

    def test_legacy_list_table_metadata_remains_supported_during_conversion(self) -> None:
        lines = ["* - Type", "  - Integer", "* - Default", "  - `2`"]
        self.assertEqual(checker.metadata_from_section(lines), ("Integer", "`2`"))

    def test_numeric_and_boolean_fortran_spellings_compare_semantically(self) -> None:
        self.assertEqual(checker.normalized_default("`1.0D0`"), checker.normalized_default("1"))
        self.assertEqual(checker.normalized_default("`.true.`"), checker.normalized_default("true"))

    def test_changed_documented_default_is_an_error(self) -> None:
        code = checker.CodeParameter("example", None, "Integer", "2")
        document = checker.DocumentedParameter("example", "Example", 5, "Integer", "`3`")
        errors, warnings = checker.check(
            [code], [document], [], [], checker.DEFAULT_DOCUMENTATION, False
        )
        self.assertEqual(len(errors), 1)
        self.assertIn("default", errors[0])
        self.assertFalse(warnings)

    def test_console_inventory_uses_aligned_columns_without_line_numbers(self) -> None:
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            checker.print_parameter_table(self.parameters[:2])
        table = output.getvalue()
        self.assertIn("Name", table)
        self.assertIn("Type", table)
        self.assertIn("Default", table)
        self.assertNotIn("\t", table)
        self.assertNotIn("Parser line", table)


if __name__ == "__main__":
    unittest.main()
