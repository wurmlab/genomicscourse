Using SSH on windows is slightly different, as windows does not come with software to do this by default. You'll need to install a client before you can connect.

## Installing MobaXTerm
There are a number of SSH clients you can use, but for the purposes of this practical, we'd recommend MobaXTerm.


1) Download from
https://mobaxterm.mobatek.net/download-home-edition.html - select installer version
2) Unzip the folder you downloaded
3) Run the file ending in .msi - you'll get a windows prompt to install this, which it is safe to agree to. Install MobaXTerm in the default location, and you're then ready to log in

## Logging in through MobaXTerm
1) Open MobaXTerm from the start menu (if you can't find it, use the search bar)


2) Click on Session in the upper left of MobaXTerm ![](images/mobasession.PNG)


3) Select SSH in the window that opens


4) In the next window, enter in the name of the server, into the Server Name box. If you've been given the server name in the form yourname@yourname.genomicscourse.com, you just need the part after the @ (so, yourname.genomicscourse.com)


5) Tick the "use default username" box, and add your username (included in the email that we sent you with log in details)


6) Tick OK, at the bottom of the box


7) You may see a message the first time you log in, saying that "the authenticity of the host can not be established", and asking if you wish to continue connecting. This message is normal when you log in to a new server. Agreeing will store the server's fingerprint, and the message will not appear again.
![](images/warning_message.png)


8) Type yes and hit return, if the message appears


9) You will then be asked for a password. Type in the password that you've been emailed, and hit return. Nothing will appear while you type. If you know you've made a mistake, you can hold down backspace, and retype the password.
![](images/linux_login_done.png)


10) Congratulations! You're now logged in!
