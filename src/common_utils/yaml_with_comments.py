#!/usr/bin/env python3
"""
YAML utilities with recursive support for commented-out values in nested structures.
"""

import os
import sys
from pathlib import Path

import ruamel.yaml
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.tokens import CommentToken
from ruamel.yaml.error import CommentMark


YamlErrorCode = 99


class CommentedData:
	"""Wrapper class to mark data that should be rendered as a comment."""
	def __init__(self, value, commented: bool = False):
		self.value = value
		self.commented = commented

	def __repr__(self):
		return f"CommentedData(value={self.value!r}, commented={self.commented})"


def format_yaml_value(value) -> str:
	"""Format a value for YAML comment representation."""
	if value is None:
		return "null"
	elif isinstance(value, (CommentedMap, dict)):
		items = []
		for k, v in value.items():
			items.append(f"{k}: {format_yaml_value(v)}")
		return "{" + ", ".join(items) + "}"
	elif isinstance(value, (CommentedSeq, list)):
		return str(list(value))
	return str(value)


def process_value(value):
	"""
	Recursively process a value, unwrapping CommentedData and converting
	nested structures to CommentedMap/CommentedSeq.

	Returns: (processed_value, is_commented)
	"""
	if isinstance(value, CommentedData):
		processed, _ = process_value(value.value)
		return (processed, value.commented)

	elif isinstance(value, dict):
		result = CommentedMap()
		active_items = []
		commented_items = []

		# First pass: process all items and categorize
		for k, v in value.items():
			processed_v, is_commented = process_value(v)
			if is_commented:
				commented_items.append((k, processed_v))
			else:
				active_items.append((k, processed_v))

		# Add active items to result
		for k, v in active_items:
			result[k] = v

		# Attach commented items as comments to last active item
		if commented_items and active_items:
			last_key = active_items[-1][0]
			comment_lines = []
			for k, v in commented_items:
				formatted = format_yaml_value(v)
				comment_lines.append(f"#  {k}: {formatted}")

			# Start with newline to force new line
			full_comment = "\n" + "\n".join(comment_lines) + "\n"
			mark = CommentMark(32)
			ct = CommentToken(full_comment, mark)
			result.ca.items[last_key] = [None, None, ct, None]

		elif commented_items and not active_items:
			# All items are commented - add as start comment
			comment_lines = []
			for k, v in commented_items:
				formatted = format_yaml_value(v)
				comment_lines.append(f"#  {k}: {formatted}")
			result.yaml_set_start_comment("\n".join(comment_lines) + "\n")

		return (result, False)

	elif isinstance(value, list):
		result = CommentedSeq()
		for item in value:
			processed, _ = process_value(item)
			result.append(processed)
		return (result, False)

	else:
		return (value, False)


def dump_data(fname: str, data: dict) -> None:
	"""
	Dump data to YAML file with support for commented-out values at any nesting level.
	"""
	path = Path(fname)
	yaml = ruamel.yaml.YAML()
	yaml.indent(mapping=2, sequence=4, offset=2)
	yaml.preserve_quotes = True

	# Process the entire data structure recursively
	processed_data, _ = process_value(data)

	try:
		with path.open('w') as fyaml:
			yaml.dump(processed_data, fyaml)
	except OSError as e:
		print(f"Cannot open/write file: {path} - {e}", file=sys.stderr)
		sys.exit(YamlErrorCode)
	except Exception as exc:
		print(f"YAML dump error: {exc}", file=sys.stderr)
		sys.exit(YamlErrorCode)


def load_data(fname: str, dbg_print: bool = False) -> dict:
	"""Load YAML data from file."""
	path = Path(fname)
	yaml = ruamel.yaml.YAML()

	try:
		with path.open('r') as fyaml:
			data = yaml.load(fyaml)
			if dbg_print:
				print(f"Loaded from {fname}:")
				print(data)
			return data if data else {}
	except OSError as e:
		print(f"Cannot open/read file: {path} - {e}", file=sys.stderr)
		sys.exit(YamlErrorCode)
	except ruamel.yaml.YAMLError as exc:
		print(f"YAML load error: {exc}", file=sys.stderr)
		sys.exit(YamlErrorCode)


# Example usage left in jump_network_builder package; the shared helper
# module intentionally does not include the execution example.


