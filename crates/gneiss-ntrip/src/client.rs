use tokio::sync::mpsc;
use url::Url;
use tracing::{info, error};
use bytes::Bytes;
use tokio::net::TcpStream;
use tokio::io::{AsyncWriteExt, AsyncReadExt};
use base64::Engine;

#[derive(Debug, Clone)]
pub struct NtripConfig {
    pub server_url: String,
    pub mountpoint: String,
    pub username: Option<String>,
    pub password: Option<String>,
}

pub struct NtripClient {
    config: NtripConfig,
}

impl NtripClient {
    pub fn new(config: NtripConfig) -> Self {
        Self { config }
    }

    /// Connects to the NTRIP caster and returns a channel receiver yielding raw bytes.
    pub async fn connect(&self) -> Result<mpsc::Receiver<Bytes>, String> {
        let url_str = format!("{}/{}", self.config.server_url.trim_end_matches('/'), self.config.mountpoint);
        let url = Url::parse(&url_str).map_err(|_| "Invalid NTRIP URL")?;

        let host = url.host_str().ok_or("Missing host in URL")?;
        let port = url.port().unwrap_or(2101);
        
        info!("Connecting to NTRIP caster: {}:{}", host, port);
        
        let mut stream = TcpStream::connect((host, port)).await.map_err(|e| e.to_string())?;

        let mut request = format!(
            "GET /{} HTTP/1.0\r\n\
             Host: {}:{}\r\n\
             User-Agent: NTRIP gneiss-client/0.1\r\n\
             Accept: */*\r\n",
            self.config.mountpoint, host, port
        );

        if let (Some(user), Some(pass)) = (&self.config.username, &self.config.password) {
            let auth = format!("{}:{}", user, pass);
            let encoded = base64::engine::general_purpose::STANDARD.encode(auth);
            request.push_str(&format!("Authorization: Basic {}\r\n", encoded));
        }
        
        request.push_str("\r\n");

        stream.write_all(request.as_bytes()).await.map_err(|e| e.to_string())?;

        // Read response header until \r\n\r\n
        let mut header_buf = Vec::new();
        let mut chunk = [0u8; 1];
        
        loop {
            let n = stream.read(&mut chunk).await.map_err(|e| e.to_string())?;
            if n == 0 {
                return Err("Connection closed by server before headers received".into());
            }
            header_buf.push(chunk[0]);
            
            if header_buf.ends_with(b"\r\n\r\n") {
                break;
            }
            if header_buf.len() > 4096 {
                return Err("Header too long".into());
            }
        }

        let header_str = String::from_utf8_lossy(&header_buf);
        if !header_str.contains("200 OK") {
            error!("NTRIP Connection failed. Response:\n{}", header_str);
            return Err("NTRIP caster rejected the connection".into());
        }

        info!("Connected successfully. Starting data stream.");

        let (tx, rx) = mpsc::channel(100);
        
        tokio::spawn(async move {
            let mut buf = [0u8; 4096];
            loop {
                match stream.read(&mut buf).await {
                    Ok(n) if n > 0 => {
                        let bytes = Bytes::copy_from_slice(&buf[0..n]);
                        if tx.send(bytes).await.is_err() {
                            break;
                        }
                    }
                    Ok(_) => {
                        info!("NTRIP stream closed cleanly by server.");
                        break;
                    }
                    Err(e) => {
                        error!("NTRIP stream read error: {}", e);
                        break;
                    }
                }
            }
        });

        Ok(rx)
    }
}
