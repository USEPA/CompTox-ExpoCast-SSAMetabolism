# -*- coding: utf-8 -*-


import pandas as pd
import pymysql
import logging
import sshtunnel
from sshtunnel import SSHTunnelForwarder

ssh_host = #USER-DEFINED
ssh_username = #USER-DEFINED
ssh_password = #USER-DEFINED

db_username = #USER-DEFINED
db_password = #USER-DEFINED

database_name = 'dev_nta_predictions_fy21'
localhost = '127.0.0.1'

def open_ssh_tunnel(verbose = False):
    if verbose:
        sshtunnel.DEFAULT_LOGLEVEL = logging.DEBUG
    
    global tunnel
    tunnel = SSHTunnelForwarder(
        (ssh_host, 22),
        ssh_username = ssh_username,
        ssh_password = ssh_password,
        remote_bind_address = ('127.0.0.1', 3306)
        )
    
    tunnel.start()
    
    return tunnel
    
def mysql_connect():
    global connection
    
    return pymysql.connect(
        host = '127.0.0.1',
        user = db_username,
        passwd = db_password,
        db = database_name,
        port = tunnel.local_bind_port)

    return connection

def run_query(sql):
    return pd.read_sql_query(sql, connection)

def mysql_disconnect():
    connection.close()
    
def close_ssh_tunnel():
    tunnel.close

