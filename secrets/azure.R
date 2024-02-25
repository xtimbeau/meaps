azure_url ="https://pourpin.blob.core.windows.net/trajets"
azure_jeton = "sp=racwdlm&st=2023-08-15T15:45:26Z&se=2024-08-15T23:45:26Z&spr=https&sv=2022-11-02&sr=c&sig=HtFyZD%2FJU%2FYqCp4WOZbBRWLIyo4vbuvkfk0g%2BIOUGR0%3D"

board <- pins::board_azure(
  AzureStor::storage_container(
    azure_url,
    sas = azure_jeton))

bd_hash <- function(obj) pins::pin_meta(board, obj)$pin_hash
bd_read <- function(obj) pins::pin_read(board, obj)

bd_write <- function(obj, name=NULL, title=NULL, description=NULL, metadata = NULL, tags=NULL, versioned=NULL) {
  if(is.null(name))
    name <- rlang::as_name(rlang::enquo(obj))
  pins::pin_write(board = board,
                  x = obj,
                  name = name, 
                  title = title,
                  description = description,
                  metadata = metadata,
                  tags = tags,
                  versioned = versioned,
                  type = "qs")
}