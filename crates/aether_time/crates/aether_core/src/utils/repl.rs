#![doc = include_str!("../../docs/repl/repl.md")]
// Some help text
const REPL_HELP_TEXT: &str = include_str!("../../docs/repl/help.txt");

use crate::parse_and_eval;
use rustyline::DefaultEditor;

pub struct Repl {}

impl Repl {
    pub fn run() {
        // This looks messed up but it is not.
        println!(
            r" _____ ______   ________  _________  ________  ___     ___    ___ 
|\   _ \  _   \|\   __  \|\___   ___\\   __  \|\  \   |\  \  /  /|
\ \  \\\__\ \  \ \  \|\  \|___ \  \_\ \  \|\  \ \  \  \ \  \/  / /
 \ \  \\|__| \  \ \   __  \   \ \  \ \ \   _  _\ \  \  \ \    / / 
  \ \  \    \ \  \ \  \ \  \   \ \  \ \ \  \\  \\ \  \  /     \/  
   \ \__\    \ \__\ \__\ \__\   \ \__\ \ \__\\ _\\ \__\/  /\   \  
    \|__|     \|__|\|__|\|__|    \|__|  \|__|\|__|\|__/__/ /\ __\ 
                                                      |__|/ \|__| 
 ________  ________           ___       ________  ________        
|\   __  \|\   ____\         |\  \     |\   __  \|\   __  \       
\ \  \|\  \ \  \___|_        \ \  \    \ \  \|\  \ \  \|\ /_      
 \ \   _  _\ \_____  \        \ \  \    \ \   __  \ \   __  \     
  \ \  \\  \\|____|\  \        \ \  \____\ \  \ \  \ \  \|\  \    
   \ \__\\ _\ ____\_\  \        \ \_______\ \__\ \__\ \_______\   
    \|__|\|__|\_________\        \|_______|\|__|\|__|\|_______|   
             \|_________|                                         "
        );

        println!("Vector Example: [1 2 3] + [4 5 6], or [1 2 3 4 5] + [5 4 3 2 1], or [1 2 3] * 2");
        println!("Matrix Example: [1 2 3; 4 5 6; 7 8 9] + [9 8 7; 6 5 4; 3 2 1]");
        println!("Type 'exit' to quit.\n");

        let mut rl = DefaultEditor::new().unwrap();

        loop {
            let readline = rl.readline(">> ");
            match readline {
                Ok(line) => {
                    let input = line.trim();
                    if input.eq_ignore_ascii_case("exit") || input.eq_ignore_ascii_case("quit") {
                        println!("Exited.");
                        break;
                    }
                    if input.eq_ignore_ascii_case("help") || input == "?" {
                        println!("{}", REPL_HELP_TEXT);
                        continue;
                    }
                    if input.is_empty() {
                        continue;
                    }
                    match parse_and_eval(input) {
                        Ok(result) => println!("{}", result),
                        Err(e) => println!("Error: {}", e),
                    }
                }
                Err(_) => break,
            }
        }
    }
}
