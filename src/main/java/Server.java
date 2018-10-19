import javax.websocket.*;
import javax.websocket.server.ServerEndpoint;

@ServerEndpoint("/")
public class Server {
    public Server(){}

    @OnOpen
    public void onOpen(Session session) {}

    @OnMessage
    public void onMessage(String message, Session session) {

    }

    @OnError
    public void onError(Throwable e){e.printStackTrace();}

    @OnClose
    public void onClose(Session session) {}
}
