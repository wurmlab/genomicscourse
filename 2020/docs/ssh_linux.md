# SSH on Linux


**Note:** Instructions are provided for Ubuntu, but if you're running a different flavour of linux, the only thing that should be different is how to open the terminal. From step 3 onwards, everything should be identical.


1) On Ubuntu, click on the icon in the top left of the screen


2) Type "Terminal", and open the terminal application

![](images/linux_find_terminal.png)


3) Into the terminal that appears, type  
   ``ssh username@servername``  
   replacing username and servername with the username and server name that you've been emailed


 ![](images/ssh_linux_login.png)


4) You will see a message the first time you log in, saying that "the authenticity of the host can not be established", and asking if you wish to continue connecting. This message is normal when you log in to a new server. Agreeing will store the server's fingerprint, and the message will not appear again.  

![](images/warning_message.png)

5) Type yes and hit return


6) You will then be asked for a password. Type in the password that you've been emailed, and hit return. Nothing will appear while you type. If you know you've made a mistake, you can hold down backspace, and retype the password.


7) Congratulations! You're now logged in! Note the change of username and server on the left of your command prompt - this tells you which server you are logged into


![](images/linux_login_done.png)