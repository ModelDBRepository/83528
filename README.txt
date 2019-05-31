Installation
============

The simulator has been tested on a couple of UNIX/Linux systems
and on Windows with cygwin. It requires something equivalent to
g++ for compilation,
lex and yacc to generate the parser of the language file and
Tcl/Tk for the result viewer look.tcl.

Get the file coco06.zip.
Then unzip, which creates the subdirectory coco06/ with the simulator in it.
Type ./make in coco06/, which creates the executable cococo
(and also prae.exe as a pre-processor).


Program call
============

Start the program in coco04/ with the following call:

./cococo -file v/mylanguagefile

Eg. mylanguagefile = HopfieldContinuous. Make sure a directory for the
activations and weights exists (in the examples, this is: /tmp/coco/).
This will start the program with it reading the language file v/mylanguagefile.
The program will allocate memory for the network you defined,
it will run the activation and training algorithm over and over again,
as you have described in your language file.

Start the ``GUI'', an image file viewer, like:

./look.tcl a w 0 1

Make sure to start this after the first weight- and activation-files have been
written out.
The arguments "a" means display activation files, "w" means display weights.
The numbers specify the areas from which to collect these files.
In the window, click the left mouse button or Return to reload the files,
right mouse button to quit.


Further reading
===============

See the file README.ps which is the same as README.pdf or docu/handbook.pdf


Background
==========

With this neural net simulator, you edit your own network and training
algorithm using a built-in ``batch''-language.
In this language, you specify every detail of your algorithm.
You can extend the language by adding and compiling to the simulator
your own functions written in C.
A variety of networks with connectionist neurons can be programmed,
in particular if neuronal updates are local.

This simulator was developed to bridge the gap between pure C-code that
becomes messy over time and simulators (such as SNNS) which restrict the user
too much. It is meant as a scheme to organize new C-code that piles up during
development. New versions are not compatible with older ones, sorry.

The need for it arose from testing new neural network learning algorithms with
different architectures, in order to explain cortico-cortical connections
(hence the name). See the publication: Emergence of modularity within one
sheet of intrinsically active stochastic neurons. C. Weber and K. Obermayer.
Proc. ICONIP, 732-737 (2000).


Content
=======

Currently, the following files are supplied, addressing the following papers:

- v/HopfieldContinuous

  This addresses no paper in particular, but such a network is used in:
  "Self-Organization of Orientation Maps, Lateral Connections, and Dynamic
   Receptive Fields in the Primary Visual Cortex"
   C. Weber, In: Proc. ICANN, 1147-52, (2001).

  Such a network also constitutes the lateral, predictive weights in
  the following:

- v/CortexDocking
  v/CortexDocking_analyze
  v/flowfield.tcl

  The first produces, the latter analyze, the toy example weights of:
  "A hybrid generative and predictive model of the motor cortex"
   Cornelius Weber, Stefan Wermter and Mark Elshaw
   Neural Networks, 19 (4), 339-353 (2006).

- v/SigmaPi_1dim

   This creates Fig.6a) of:
   "A Self-Organizing Map of Sigma-Pi Units",
    C. Weber and S. Wermter, Neurocomputing (2006/7, in press).

- v/Vert_Horiz_Saccades

   This creates Fig.5 of:
   "A possible representation of reward in the learning of saccades"
   C. Weber and J. Triesch, Proc. Epigenetic Robotics, pp. 153-60 (2006).


Correspondence
==============

Dr. Cornelius Weber, Room 0.318, Frankfurt Institute for Advanced Studies,
Johann Wolfgang Goethe University, Max-von-Laue Str. 1, 60438 Frankfurt am
Main, Germany. Tel: +49 69 798 47536, Fax: +49 69 798 47611.
WWW: http://fias.uni-frankfurt.de/~cweber/
Email: c.weber@fias.uni-frankfurt.de
