



count = 0

n = 1

global opsin

opsin = ""

with open("anc_opsins.txt","r") as file:
    opsin = ""

    lines = file.readlines()
    i = int(len(lines))
    while count <= 12:
    
        opn = n
        content = str(lines[opn-1].strip())
        species = content.split('(')    
        name = str(species[0])

        while opn <= i:
        
            content = str(lines[opn]).strip()
            content = content[6:]
            opn += 26
            opsin = opsin +content
       
        print("Ancient Opsin " + name + ": " + opsin.strip())

        n += 2
        count += 1
        opsin = ""

    file.close()