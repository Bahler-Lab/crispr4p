#!/usr/local/bin/python2.7


# Import modules for CGI handling
import cgi, cgitb
import os
import crispr4p
from flup.server.fcgi import WSGIServer

print("Content-type: text/html\n\n")

def load_bahler_template():
    src_path = os.path.dirname(__file__)
    try:
        template_file = open(src_path + '/bahler_template.html')
        print template_file.read()
    except IOError as err:
        print err
#
#
#form = cgi.FieldStorage()
#print """
#<form action="app" method="post">
#First Name: <input type="text" name="first_name"><br />
#Last Name: <input type="text" name="last_name" />
#
#<input type="submit" value="Submit" />
#</form> """
#

def app(env, start_response):
    start_response('200 OK', [('Content-Type', 'text/plain')])
    load_bahler_template()
    form = cgi.FieldStorage()
    yield ' '
    generate_form()
    if form.has_key("action") and form.has_key("name"):
        if (form["action"].value == "display"):
            display_data(form["name"].value, form["age"].value)
            yield 'xx'
            #os.system("/var/www/bgi-bin/cris4p/crispr4p.py -cr I -co 112000...115000")




# Define function to generate HTML form.
def generate_form():
    #print "<HTML>\n"
    #print "<HEAD>\n"
    #print "\t<TITLE>CRIS4P</TITLE>\n"
    #print "</HEAD>\n"
    #print "<BODY BGCOLOR = white>\n"
    print "<H3>Please, enter gene information.</H3>\n"
    print "<TABLE BORDER = 0>\n"
    print "<FORM METHOD = post ACTION = \"cris4p_web.py\">\n"
    print "\t<TR><TH>Name:</TH><TD><INPUT type = text name = \"name\"></TD><TR>\n"
    print "\t<TR><TH>Chromosome:</TH><TD><INPUT type = text name = \"chromosome\"></TD></TR>\n"
    print "\t<TR><TH>Coordinate:</TH><TD><INPUT type = text name = \"coor_lower\"></TD>"
    print "<TD><INPUT type = text name = \"coor_upper\"></TD>\n"
    print "</TABLE>\n"
    print "<INPUT TYPE = hidden NAME = \"action\" VALUE = \"display\">\n"
    print "<INPUT TYPE = submit VALUE = \"Enter\">\n"
    print "</FORM>\n\n"
    #print "</BODY>\n"
    #print "</HTML>\n"

    # Define function display data.
def display_data(name):
    print "<HTML>\n"
    print "<HEAD>\n"
    print "\t<TITLE>Info Form</TITLE>\n"
    print "</HEAD>\n"
    print "<BODY BGCOLOR = white>\n"
    print name, ", you are", age, "years old."
    print "</BODY>\n"
    print "</HTML>\n"

    # Define main function.
def main():
    load_bahler_template()
    form = cgi.FieldStorage()
    generate_form()
    if form.has_key("action") and form.has_key("name"):
        if (form["action"].value == "display"):
            print form["name"].value
            #display_data(form["name"].value, form["age"].value)

    # Call main function.
main()

#WSGIServer(app).run()

