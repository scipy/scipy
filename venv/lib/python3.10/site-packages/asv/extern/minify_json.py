'''
Created on 20/01/2011

v0.2 (C) Gerald Storer
MIT License

Based on JSON.minify.js:
https://github.com/getify/JSON.minify

Contributers:
 - Pradyun S. Gedam (conditions and variable names changed)
'''


import re

def json_minify(string, strip_space=True):
    tokenizer = re.compile(r'"|(/\*)|(\*/)|(//)|\n|\r')
    in_string = False
    in_multi = False
    in_single = False

    new_str = []
    index = 0

    for match in re.finditer(tokenizer, string):

        if not (in_multi or in_single):
            tmp = string[index:match.start()]
            if not in_string and strip_space:
                # replace white space as defined in standard
                tmp = re.sub('[ \t\n\r]+', '', tmp)
            new_str.append(tmp)

        index = match.end()
        val = match.group()

        if val == '"' and not (in_multi or in_single):
            escaped = re.search(r'(\\)*$', string[:match.start()])

            # start of string or unescaped quote character to end string
            if not in_string or (escaped is None or len(escaped.group()) % 2 == 0):
                in_string = not in_string
            index -= 1 # include " character in next catch
        elif not (in_string or in_multi or in_single):
            if val == '/*':
                in_multi = True
            elif val == '//':
                in_single = True
        elif val == '*/' and in_multi and not (in_string or in_single):
            in_multi = False
        elif val in '\r\n' and not (in_multi or in_string) and in_single:
            in_single = False
        elif not ((in_multi or in_single) or (val in ' \r\n\t' and strip_space)):
            new_str.append(val)

    new_str.append(string[index:])
    return ''.join(new_str)


if __name__ == '__main__':
    # Python 2.6+ needed to run tests
    import json
    import textwrap
    import unittest

    class JsonMinifyTestCase(unittest.TestCase):
        """Tests for json_minify"""
        def template(self, in_string, expected):
            in_dict = json.loads(json_minify(in_string))
            expected_dict = json.loads(expected)
            self.assertEqual(in_dict, expected_dict)

        def test_1(self):
            self.template(textwrap.dedent('''
                // this is a JSON file with comments
                {
                    "foo": "bar",    // this is cool
                    "bar": [
                        "baz", "bum"
                    ],
                /* the rest of this document is just fluff
                   in case you are interested. */
                    "something": 10,
                    "else": 20
                }

                /* NOTE: You can easily strip the whitespace and comments
                   from such a file with the JSON.minify() project hosted
                   here on github at http://github.com/getify/JSON.minify
                */'''),
                '{"foo":"bar","bar":["baz","bum"],"something":10,"else":20}'
            )

        def test_2(self):
            self.template(textwrap.dedent('''
                {"/*":"*/","//":"",/*"//"*/"/*/"://
                "//"}'''),
                '{"/*":"*/","//":"","/*/":"//"}'
            )

        def test_3(self):
            self.template(textwrap.dedent(r'''
                /*
                this is a
                multi line comment */{

                "foo"
                :
                    "bar/*"// something
                    ,    "b\"az":/*
                something else */"blah"

                }
                '''),
                r'{"foo":"bar/*","b\"az":"blah"}'
            )

        def test_4(self):
            self.template(textwrap.dedent(r'''
                {"foo": "ba\"r//", "bar\\": "b\\\"a/*z",
                "baz\\\\": /* yay */ "fo\\\\\"*/o"
                }
                '''),
                r'{"foo":"ba\"r//","bar\\":"b\\\"a/*z","baz\\\\":"fo\\\\\"*/o"}'
            )

    unittest.main()
