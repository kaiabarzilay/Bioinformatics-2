#!/usr/bin/env python
# basic_classes.py

class Circle():
	''' A class called Circle that uses radius and color. '''
	def __init__(self, color, radius):
		self.color = color
		self.radius = radius

	def diameter(self):
		''' Find the diameter of the circle '''
		return 2 * self.radius

	def circumference(self):
		''' Find the circumference of the circle '''
		return 2 * 3.14 * self.radius

	def isRed(self):
		''' Return True is circle is red '''
		if self.color == "red":
			return color


circle_1 = Circle("blue", 10)
print circle_1.circumference()

class GraduateStudent():
	''' A class called GraduateStudent that uses names, year and major '''
	def __init__(self, first_name, last_name, year, major):
		self.first_name = first_name
		self.last_name = last_name
		self.year = year
		self.major = major

	def year_matriculated(self):
		''' return the matriculation year '''
		return 2020 - self.year

matriculation_year = GraduateStudent("Nancy", "Smith", 2017, "Buisness")
print matriculation_year.year_matriculated()	
