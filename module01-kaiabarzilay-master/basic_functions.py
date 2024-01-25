#!/usr/bin/env python
# basic_functions.py

def multiply(a, b):
	"""
	return the product of two numbers
	"""
	return a * b

# Test of integers
x = multiply(5, 5)
print(x)

def hello_name(string=""):
	"""
	Return "Hello, string!" or "Hello, you!"
	"""
	if string=="":
		return "Hello, you!"
	else:
		return "Hello," + string + "!"

blank = hello_name()
print(blank)
# Test of a name
name = hello_name("Kaia")
print(name)


def less_than_ten(sequence):
	""" 
	Return numbers in list < 10
	"""
	numbers = []
	for x in sequence:
		if x < 10:
			numbers.append(x)
	return numbers

# Test of output
numbers = less_than_ten([1, 5, 81, 10, 8, 2, 102])
print(numbers)
