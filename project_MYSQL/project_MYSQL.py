import mysql.connector



def create_connection(
        inp_host_name,
        inp_port,
        inp_db_name,
        inp_user_name,
        inp_password):
    '''
        Description:
            Creates a MYSQL connection
        
        Inputs:
            inp_host_name: str: The host
            inp_port: int: The port
            inp_db_name: str: The DB name
            inp_user_name: str: The user name
            inp_password: str: The user password
        
        Outputs:
            out_connection: : The MYSQL connection

    '''
    print ("aaaab")
    mydb = mysql.connector.connect(
        host = inp_host_name,
        port = inp_port,
        user = inp_user_name,
        password = inp_password,
        database = inp_db_name)



if __name__ == "__main__":
    host_name = "localhost"
    port = 3306
    db_name = "sakila"
    user_name = "local_aee"
    password = "KargaCrow01"
    conn = create_connection(
        host_name,
        port,
        db_name,
        user_name,
        password)
    print (conn)


