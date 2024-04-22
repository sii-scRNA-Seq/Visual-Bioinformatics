import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { DomSanitizer } from '@angular/platform-browser';
import { fakeAsync, tick, TestBed } from '@angular/core/testing';
import { first, firstValueFrom } from 'rxjs';
import { HttpClientTestingModule, HttpTestingController } from '@angular/common/http/testing';
import { MatSnackBar, MatSnackBarModule } from '@angular/material/snack-bar';

import { Block } from './block.interface';
import { MockUserIdService } from './mock-user-id.service';
import { OutputService } from './output.service';
import { UserIdService } from './user-id.service';

import { createServer, Server } from "http";
import { AddressInfo } from "net";
import { io as ioc, Socket } from "socket.io-client";
import { BackendSocketClient } from './backend-socket.client';
import { MockBackendSocketClient } from './mock-backend-socket.client';

describe('OutputService', () => {
  let service: OutputService;
  let snackBar: MatSnackBar;
  let sanitizer: DomSanitizer;
  let io: Server, serverSocket: Socket;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        BrowserAnimationsModule,
        HttpClientTestingModule,
        MatSnackBarModule,
      ], 
      providers: [
        { provide: UserIdService, useClass: MockUserIdService },
        { provide: BackendSocketClient, useClass: MockBackendSocketClient }
      ],
    });
    service = TestBed.inject(OutputService);
    snackBar = TestBed.inject(MatSnackBar);
    sanitizer = TestBed.inject(DomSanitizer);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('outputs', () => {
    it('should initially be empty', async () => {
      const outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(0);
    });
  });

  describe('resetOutputs', () => {
    it('should set the outputs list to be empty', async () => {
      let outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(0);
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      const blocks: Block[] = [];
      service.executeBlocks(blocks);
      outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(1);
      service.resetOutputs();
      outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(0);
    });
  });

  describe('executeBlocks', () => {
    it('should open snack bar when userId is null', () => {
      const spy = spyOn(snackBar, 'open');
      const blocks: Block[] = [];
      service.executeBlocks(blocks);
      expect(spy).toHaveBeenCalledOnceWith('No User ID, please refresh the page and try again', 'Close', { duration: 5000 });
    });

    it('should call backendSocketClient.sendRequest() with the appropriate message', () => {
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      const backendSocketClient: BackendSocketClient = TestBed.inject(BackendSocketClient);
      const spy = spyOn(backendSocketClient, 'sendRequest');
      const block: Block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {key: 'foo', text: 'bar', value: 42},
        ],
      };
      const blocks: Block[] = [];
      blocks.push(block);
      service.executeBlocks(blocks);
      const message = {
        user_id: 'mock_user_id',
        blocks: [{ block_id: 'loaddata', foo: 42}],
      };
      expect(spy).toHaveBeenCalledOnceWith(message);
    });
  });

  describe('backendSocketClient subscription', () => {
    it('should add text to outputs array when it receives a valid text response', async () => {
      let outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
      const backendSocketClient: BackendSocketClient = TestBed.inject(BackendSocketClient);
      backendSocketClient.sendRequest("text");
      outputs = await firstValueFrom(service.outputs);
      expect(outputs.length).toBe(1);
      expect(outputs[0].text).toBe('Dummy text');
    });

    it('should add sanitised SafeUrl image and alt text to outputs array when it receives a valid image response', async () => {
      let outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
      const backendSocketClient: BackendSocketClient = TestBed.inject(BackendSocketClient);
      const spy = spyOn(sanitizer, 'bypassSecurityTrustUrl');
      backendSocketClient.sendRequest("image");
      expect(spy).toHaveBeenCalledTimes(1);
      outputs = await firstValueFrom(service.outputs);
      expect(outputs.length).toBe(1);
      const expectedValue = sanitizer.bypassSecurityTrustUrl('data:image/png;base64,' + 'Image string');
      expect(outputs[0].img).toBe(expectedValue);
      expect(outputs[0].alttext).toBe('Alt text');
    });

    it('should not change the outputs array when it receives an end_connection response', async () => {
      let outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
      const backendSocketClient: BackendSocketClient = TestBed.inject(BackendSocketClient);
      backendSocketClient.sendRequest("end_connection");
      outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
    });

    it('should open snack bar with given message when response is a known error', async() => {
      let outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
      const backendSocketClient: BackendSocketClient = TestBed.inject(BackendSocketClient);
      const spy = spyOn(snackBar, 'open');
      backendSocketClient.sendRequest("error");
      expect(spy).toHaveBeenCalledOnceWith('error', 'Close', { duration: 5000 });
      outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
    });

    it('should open snack bar with generic message for invalid responses', async() => {
      let outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
      const backendSocketClient: BackendSocketClient = TestBed.inject(BackendSocketClient);
      const spy = spyOn(snackBar, 'open');
      backendSocketClient.sendRequest("invalid response");
      expect(spy).toHaveBeenCalledOnceWith('Received a bad response, please refresh the page and try again', 'Close', { duration: 5000 });
      outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
    });
  });

  describe('executingBlocks', () => {
    it('should work', async () => {
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      expect(false).toBe(true);
    });
  });
});
