/** 
    @file hydrodoc.h
    @brief Main page of the documentation
 
    @par CREATION
    @author Stefan Schippers
    
    @par VERSION
    $Id: hydrodoc.h 372 2016-02-04 18:43:22Z iamp $

 */ 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @mainpage
 *
 * @author 
 *  Stefan Schippers<br />
 *  Atomic and Molecular Physics<br />
 *  I. Physics Institute<br />
 *  Justus-Liebig-University Giessen<br />
 *  Leihgesterner Weg 217<br />
 *  35392 Giessen, Germany<br />
 *<br />
 *  Stefan.Schippers@physik.uni-giessen.de
 *
 * @include doc/svnversion.h
******************************************************************************************************************
 * @section Introduction
 *
 * The program hydrocal is a collection of subroutines performing hydrogenic atomic structure calculations such as
 * bound-bound and bound-free radiative transion rates which are used in, e.g., the calaculation of cross sections
 * and rate coefficients for radiative recombination. In addition convolutions can be carried out, e.g., for the 
 * calculation of rate coefficients from cross sections. The program suite has evolved over the years. The development was
 * driven mainly by data evaluation needs of storage-ring electron-ion collision experiments.
 * Some of the calculations are explained in the following references
 * -# S. Schippers et al. <a href="http://dx.doi.org/10.1088/0953-4075/28/15/017">JPB 28 (1995) 3271</a>
 * -# S. Schippers et al. <a href="http://dx.doi.org/10.1086/321512">ApJ 555 (2001) 1027</a> 
 * -# S. Schippers et al. <a href="http://dx.doi.org/10.1051/0004-6361:20040380">A&A 421 (2004) 1185</a> 
 *
******************************************************************************************************************
 * @section Installation 
 *
 * All files that are required for the installation can be retrieved from the svn repository
 * http://www.strz.uni-giessen.de/repository. 
 * Documentation about svn is provided at http://svnbook.red-bean.com/index.en.html.
 * For checking out the hydrocal source files change to your home directory and use the following command:
 * @verbatim
    svn checkout http://www.strz.uni-giessen.de/repository/hydrocal hydrocal
@endverbatim
 *
 * Required Linux packages (all these packages are part of any major Linux distribution): 
 * -# <a href="https://www.gnu.org/software/make/">GNU make</a>,
 * -# <a href="https://gcc.gnu.org/onlinedocs/gcc/index.html">GNU C++ compiler</a> (g++), 
 * -# <a href="http://www.stack.nl/~dimitri/doxygen/index.html">doxygen</a> with (optionally) 
 * <a href="http://www.graphviz.org/">graphviz</a> for generating the documentation. 
 *
 * All software can be compiled and installed by using the <tt>make</tt> command (see @ref makefile for details). 
 * -# Go to the hydrocal main directory, usually <tt>~/hydrocal</tt>.
 * -# Type <tt>make</tt> to obtain a list of available options.
 * -# Type <tt>make hydrocal</tt> in order to compile the program itself.
 * -# Type <tt>make install</tt> to install the program in the ~/bin folder.
 * -# Type <tt>make install-html</tt> to generate the html documentation.
 * -# Type <tt>make install-pdf</tt> to generate a pdf file containing the entire documentation.
 * 

******************************************************************************************************************
 * @section Usage
 *
 * The program can be invoked by simply typing <tt>hydrocal</tt> on the command line. There are no command-line parameters.
 * The program prompts for all necessary input. There are several groups of subprograms which can be selected 
 * from a tree-like menu structure (just explore or see @ref hydrocal.cxx for details).
 */

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** @page makefile makefile
 *
 *  Type <tt>make</tt> on the command line to see a list of the available options
 *
 *  @verbinclude makefile
 *
 */
