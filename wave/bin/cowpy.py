#!/usr/bin/env python

# encoding: utf-8

from cowpy import cow

# Create a Cow
cheese = cow.Moose()

# Get a cowsay message by milking the cow
msg = cheese.milk("My witty message")

# do something with the message
print(msg)
