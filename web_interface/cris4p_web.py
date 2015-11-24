#!/usr/bin/python

# Import modules for CGI handling
import cgi, cgitb

# Create instance of FieldStorage
form = cgi.FieldStorage()

print "Content-type: text/html\n\n"

print ""




#print
#    start_form(),
#
#    p(qq(The help file for this program is available <A HREF='../../PPPP/help_file_deletion.shtml' TARGET='resource window'> here</A>) ),
#    p(qq(Enter the <B>name of gene</B> to delete @{[textfield('gene')]}) ),
#    p(qq(Enter desired <B>length of primer target sequence</B> @{[textfield('length', '80', 3)]}) ),
#    p(qq(i.e. the length of primer excluding plasmid specific sequence.) ),
#    p(qq(Which plasmid will you use as a <B>PCR template?</B>) ),
#    radio_group(-name    => 'plasmid',
#		-default => 'pFA6a',
#		-values  => ['pFA6a', 'KS-ura4', 'Other']),
#
#    p(qq(Five pairs of primers will be suggested to you. The primers are directly
#	 upstream of the ORF (not including the ORF sequence). The remaining pairs are positioned
#	 away from the ORF in user-defined increments.
#	 How far away (base pairs from ORF) would you like each subsequent primer set to be?
#	 Minus values acceptable, which will result in primers encroaching into the ORF.) ),
#    p(qq(<B> Primer increment</B>  @{[textfield('increment', '40', 3)]}) ),
#    p(qq(Long strings of identical bases are highlighted in<FONT COLOR ='red'> colour</FONT>.) ),
# 
#    submit(),
#    reset(),
#    end_form(),
# p(qq(</td> <td><img src="/gfx/pppp.jpg" alight=right>
#</td></tr><tr>)),
# p(qq(</tr></table>)),
#

