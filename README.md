

**VenusTools** is a set of functions designed to help working with Venus data from the Venus Express mission.

There are currently two different pseudo Python classes in **VenusTools**, which are stand-alone (they do not depend on each other): **VeRaTools** and **VMCTools**. 

There is no ```pip install``` option (yet).
Just clone this repository to a directory of your choice on your computer.
In order to use these classes, the Python session must know the paths to the directories on your machine.
Hence you need to tell Python where to look. If you use a Python startup file the following paths need to be appended to the ```sys.path``` variable:

  ```
  sys.path.append ('/SomeWhereOnYourMachine/VenusTools/VeRaTools')
  sys.path.append ('/SomeWhereOnYourMachine/VenusTools/VMCTools')
  ```

If you run Python straight from the command line, then on a Mac you need to add the following lines to .zprofile (or .bashrc or similar):

```
export PYTHONPATH="/SomeWhereOnYourMachine/VenusTools/VeRaTools:$PYTHONPATH"
export PYTHONPATH="/SomeWhereOnYourMachine/VenusTools/VMCTools:$PYTHONPATH"
```

The methods in these classes also use methods from the [GeneralTools]{https://github.com/PleaseStateTheNatureOfYourInquiry/GeneralTools},
also developed by Maarten Roos-Serote.

Please find the documentation on how to use **VenusTools** [here](https://venusresearchworkbook.readthedocs.io/en/latest/index.html).
