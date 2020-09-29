function newstr = underline2space(str)

newstr = str;
newstr(str=='_') = ' ';
