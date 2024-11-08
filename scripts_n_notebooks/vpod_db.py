
"""
vpod_db.py : routines for accessing vpod database using sqlite3

init_db(db_path, data_path=None) - opens the database file, initializing it with
data from 'data_path' if the database is empty.
"""

import sqlite3
from pathlib import Path

def open_db(db_path):
    return sqlite3.connect(db_path)

def check_db(db):
    '''check to see if the given database has the required tables'''
    try:
        #attempt to read one entry from each of the tables
        # an exception will be thrown if the table does not exist
        db.execute('SELECT * FROM scp LIMIT 1;')
        db.execute('SELECT * FROM heterologous LIMIT 1;')
        db.execute('SELECT * FROM links LIMIT 1;')
        db.execute('SELECT * FROM litsearch LIMIT 1;')
        db.execute('SELECT * FROM opsins LIMIT 1;')
        db.execute('SELECT * FROM refs LIMIT 1;')
        return True
    except:
        return False

def init_db(db_path, data_path=None):
    '''Open the database at db_path, creating it if necessary.
    If the database does not contain the needed tables, it will be populated by
    the data at 'data_path'
    '''
    db = open_db(db_path)
    if not check_db(db) and data_path is not None:
        create_schema(db)
        import_literature(db, data_path)
        import_references(db, data_path)
        import_heterologous(db, data_path)
        import_opsins(db, data_path)
    return db

def create_schema(db):
    cur = db.cursor()

    #try:
    #    cur.execute("""
    #    DROP DATABASE vizphiz_db;
    #    """)
    #    db.commit() 
    #except:
    #    "vizphiz_db does not yet exist!"
    #    pass

    #cur.execute("""CREATE DATABASE vizphiz_db;""")
    #db.commit() 

    #cur.execute("""USE vizphiz_db;""")
    #db.commit() 


    cur.execute("""
    CREATE TABLE scp
    (
    id int unsigned not null primary key,
    genus varchar(50),
    species varchar(50),
    celltype varchar(50),
    cellsubtype varchar(50),
    lamdamax decimal(9,5),
    error decimal(9,5),
    chromophore varchar(50),
    method varchar(50),
    stage varchar(50),
    refid int,
    notes varchar(1000)
    );
    """)
    db.commit() 


    cur.execute("""
    CREATE TABLE heterologous
    (
    hetid int unsigned not null primary key,
    genus varchar(50),
    species varchar(50),
    accession varchar(500),
    mutations varchar(500),
    lamdamax decimal(9,5),
    error decimal(9,5),
    cellculture varchar(50),
    purification varchar(50),
    spectrum varchar(50),
    sourcetype varchar(50),
    refid int,
    notes varchar(1000)
    );
    """)
    db.commit() 


    cur.execute("""
    CREATE TABLE links
    (
    linkid int unsigned not null primary key,
    accession varchar(500),
    maxid int,
    refid int,
    evidence varchar(1000)
    );
    """)
    db.commit() 


    cur.execute("""
    CREATE TABLE litsearch
    (
    searchid int,
    researcher varchar(50),
    month int,
    year int,
    engine varchar(500),
    keywords varchar(500)
    );
    """)
    db.commit() 


    cur.execute("""
    CREATE TABLE opsins
    (
    opsinid int unsigned not null primary key,
    genefamily varchar(50),
    genenames varchar(50),
    genus varchar(50),
    phylum varchar(25),
    class varchar(25),
    species varchar(50),
    db varchar(50),
    accession varchar(500),
    dna varchar(10000),
    aa varchar(3333),
    refid int
    );
    """)

    db.commit()

    cur.execute("""
    CREATE TABLE refs
    (
    refid int,
    doi varchar(1000),
    searchid int
    );
    """)
    db.commit()

def import_literature(db, path):
    with open(Path(path)/'litsearch.csv','r') as f:
        lines = f.readlines()
    
    count=0
    for line in lines:
        if count == 0:
            count+=1
        else:
            columns = line.split("\t")

            mycursor = db.cursor()

            sql = "INSERT INTO vizphiz_db.litsearch (searchid, researcher, month, year, engine, keywords) VALUES (%s, %s, %s,%s, %s, %s)"
            val = (columns[0], columns[1], columns[2], columns[3], columns[4], columns[5])
            print(sql)
            print(val)

            mycursor.execute(sql, val)

            db.commit()

            print(mycursor.rowcount, "record inserted.")

def import_references(db, data_path):
    with open(Path(data_path)/'references.csv','r') as f:
        lines = f.readlines()
    count=0
    for line in lines:
        if count == 0:
            count+=1
        else:
            columns = line.split("\t")

            mycursor = db.cursor()

            sql = "INSERT INTO vizphiz_db.refs (refid, doi, searchid) VALUES (%s, %s, %s)"
            val = (columns[0], columns[1], columns[2])
            print(sql)
            print(val)

            mycursor.execute(sql, val)

            db.commit()

            print(mycursor.rowcount, "record inserted.")

def import_heterologous(db, data_path):
    with open(Path(data_path)/'heterologous.csv','r',encoding='utf8') as f:
        lines = f.readlines()
    count=0
    for line in lines:
        if count == 0:
            count+=1
        else:
            columns = line.split("\t")
            print(columns)
            mycursor = db.cursor()

            sql = "INSERT INTO vizphiz_db.heterologous (hetid, genus, species, accession, mutations, lamdamax, error, cellculture, purification, spectrum, sourcetype, refid) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
            val = (columns[0], columns[1], columns[2], columns[3], columns[4], columns[5], columns[6], columns[7], columns[8], columns[9], columns[10], columns[11])
            print(sql)
            print(val)

            mycursor.execute(sql, val)

            db.commit()

            print(mycursor.rowcount, "record inserted.")

def import_opsins(db, data_path):
    with open(Path(data_path)/'opsins.csv', 'r') as f:
        lines = f.readlines()

    count=0
    for line in lines:
        if count == 0:
            count+=1
        else:
            columns = line.split("\t")

            mycursor = db.cursor()

            sql = "INSERT INTO vizphiz_db.opsins (opsinid, genefamily, genenames, phylum, class, genus, species, db, accession, dna, aa, refid) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
            val = (columns[0], columns[1], columns[2], columns[3], columns[5], columns[6], columns[7], columns[8], columns[9], columns[10], columns[11], columns [12])
            print(sql)
            print(val)

            mycursor.execute(sql, val)

            db.commit()

            print(mycursor.rowcount, "record inserted.")
