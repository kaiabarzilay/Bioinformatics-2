#!/usr/bin/env python
# hello_you.py

import sys

def hello_name(string=""):
        """
        Return "Hello, string!" or "Hello, you!"
        """
        if string=="":
                return "Hello, you!"
        else:
                return "Hello," + string + "!"

if __name__ == "__main__":
	# Check to make sure that there are at least two arguments
	arg_count = len(sys.argv) - 1
	if arg_count == 1:
		name_value = sys.argv[1]
		hello_to_print = hello_name(name_value)
		print(hello_to_print)
	else:
		no_input = hello_name()
		print(no_input)
