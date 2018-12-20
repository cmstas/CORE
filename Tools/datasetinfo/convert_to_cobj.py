lines = open("scale1fbs.txt").readlines()

content = 'const std::unordered_map<std::string, DatasetInfoFromFile::datasetInfo> default_dslist = {\n'
last_line = len(lines) - 1
for i, line in enumerate(lines[:-1]):
    if line.startswith('#'): continue
    line = line.split()
    if line[4] == "1": continue
    content += '  {{"{1}{0}", {{{2}, {3}, {4}, {5}}}}}'.format(line[0], line[1], line[2], line[3], line[4], line[5])
    if i < last_line:
        content += ',\n'

line = lines[-1]
if not line.startswith('#'):
    line = line.split()
    content += '  {{"{1}{0}", {{{2}, {3}, {4}, {5}}}}}\n'.format(line[0], line[1], line[2], line[3], line[4], line[5])
    content += '};'
else:
    print "Warning, last line not valid!!"

f = open("defaultList.icc", 'w')
f.write(content)

print "Done!"
