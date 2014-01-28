# custom test functions

is_close <- function(y) {
	function(x) {
		expectation( abs(x-y)/abs(y) <= 1e-3, "isn't close")
	}
}
